//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/filesystem/operations.hpp>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include "clj.h"

namespace LaDa
{
  namespace Models
  {
    bool Clj :: Load( const TiXmlElement& _node )
    {
      const TiXmlElement* const parent = opt::find_node( _node, "Functional", "Clj" );
      if( not parent ) return false;
      if( not Ewald::Load( *parent ) ) return false;
      if( not LennardJones :: Load( *parent ) ) return false;
      return true;
    }

    Clj :: t_Return Clj :: energy(const t_Arg& _in, t_Arg& _out) const
    {
      _out.cell.zero();
      foreach( t_Arg :: t_Atom &atom, _out.atoms )
        atom.pos = atat::rVector3d(0,0,0);

      _out.energy = Ewald::energy( _in, _out ) + LennardJones::energy( _in, _out );
      return _out.energy;
    }

    namespace details
    {
      // contains atom specific quantities.
      struct atomic_species
      {
        size_t atomic_number;
        types::t_real radius;
        types::t_real charge;
      };
      bool operator==( const bond_type &_a, const bond_type &_b )
      {
        return _a.atomic_number == _b.atomic_number;
      }
      // contains bond specific quantities.
      struct bond_type
      {
        size_t a;
        size_t b;
        types::t_real epsilon;
        types::t_real rsigma;
      };
      bool operator==( const bond_type &_a, const bond_type &_b )
      {
        return (_a.a == _b.a and  _a.b == _b.b ) or (_a.a == _b.b and  _a.b == _b.a );
      }
    }

    void read_fortran_input( Clj &_clj, boost::filesystem &_path )
    {
      __TRYBEGIN

      
      namespace bsc = boost::spirit::classic;
      namespace fs = boost::filesystem;  
      __DOASSERT( not fs::exists( _path ), "Path " << _path << " does not exits.\n" )
      __DOASSERT( not( fs::is_regular( _path ) or fs::is_symlink( _path ) ),
                  _path << " is neither a regulare file nor a system link.\n" )
      std::ifstream file( _path.string().c_str(), std::ifstream::in );
      std::string line;
      __DOASSERT( file.bad(), "Could not open file " + _path.string() + "\n" );
      
      // Number of atomic species.
      std::getline( file, line );
      __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
      const size_t nbspecies( boost::lexical_cast< size_t >(line) );

      // Atomic types
      __DOASSERT( nbspecies == 0, "Number of species set to zero in input " + _path.string() + ".\n" )
      std::vector< atomic_species > species;
      for( size_t i(0); i < nbspecies; ++i )
      {
        atomic_species atom;
        std::getline( file, line );
        __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
        __DOASSERT
        (
          not bsc::parse
          (
            line.c_str(),
          
                bsc::uint_p[ bsc::assign_a( atom.atomic_number ) ] 
             >> bsc::ureal_p[ bsc::assign_a( atom.radius ) ] 
             >> bsc::real_p[ bsc::assign_a( atom.charge ) ],
             bsc::space_p
          ).hit,
          "Could not parse atomic type.\n" 
        )
        __DOASSERT( species.end() != std::find( species.begin(), species.end(), atom ),
                    "Same specie specified twice.\n" )
        species.push_back( atom );
      }
      
      // Bond types
      __DOASSERT( nbspecies == 0, "Number of species set to zero in input " + _path.string() + ".\n" )
      std::vector< bond_type > bonds;
      for( size_t i(0); i < nbspecies * nbspecies; ++i )
      {
        bond_type bond;
        std::getline( file, line );
        __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
        __DOASSERT
        (
          not bsc::parse
          (
            line.c_str(),
          
                bsc::uint_p[ bsc::assign_a( bond.a ) ] 
             >> bsc::uint_p[ bsc::assign_a( bond.b ) ] 
             >> bsc::ureal_p[ bsc::assign_a( atom.epsilon ) ] 
             >> bsc::real_p[ bsc::assign_a( atom.rsigma ) ],
             bsc::space_p
          ).hit,
          "Could not parse atomic type.\n" 
        )
        const std::vector< bond_type > :: const_iterator i_found
        ( 
          std::find( bonds.begin(), bonds.end(), bond )
        );
        if( i_found != bonds.end() )
        {
          __DOASSERT( i_found->a == bond.a, i_found->b == bond.b,
                      "Same bond specified twice (bond " << i << ").\n") 
          __DOASSERT( std::abs( i_found->epsilon - bond.epsilon ) > 1e-12,
                      "Same bond, different epsilon (bond " << i << ").\n") 
          __DOASSERT( std::abs( i_found->rsigma - bond.rsigma ) > 1e-12,
                      "Same bond, different rsigma (bond " << i << ").\n") 
          continue;
        }
        bonds.push_back( bond );
      }

      // lj real space cutoff
      std::getline( file, line );
      __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
      LennardJones :: rcut_ = boost::lexical_types<types::t_real>( Print::StripEdges( line ) );
      __DOASSERT( rcut_ < 0e0, "Negative real space cutoff for LJ.\n" );

      // lj real space cutoff
      std::getline( file, line );
      __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
      Ewald :: cutoff_ = boost::lexical_types<types::t_real>( Print::StripEdges( line ) );
      __DOASSERT( rcut_ < 0e0, "Negative real space cutoff for Ewald.\n" );

      // PEGS.
      std::getline( file, line );
      __DOASSERT( file.eof(), "Unexpected end-of-file " + _path.string() + "\n" );
       
      __DOASSERT( boost::lexical_types<types::t_real>( Print::StripEdges( line ) ) != 1e0,
                  "No pegs factor implemented.\n" )

      // Now converts to Ewald.
      Ewald :: charges_.clear(); 
      foreach( const atomic_specie & specie, species )
      {
        __DOASSERT( Physics :: Symbol( specie.atomic_number ) == "error",
                    "Unknow atomic number " << specie.atomic_number << ", please complain.\n" )
        Ewald :: charges_[ Physics :: Symbol( specie.atomic_number ) ] = specie.charge;
      }

      // Now converts to LJ.
      LennardJones :: bonds_.clear();
      foreach( const bond_type & bond, bonds )
      {
        __DOASSERT( Physics :: Symbol( bond.a ) == "error",
                    "Unknow atomic number " << bond.a << ", please complain.\n" )
        __DOASSERT( Physics :: Symbol( bond.b ) == "error",
                    "Unknow atomic number " << bond.b << ", please complain.\n" )
        __DOASSERT( species.end() == find( species.begin(), species.end(), bond.a ), 
                    "Atomic type " << bond.a << " incomplete in input.\n" )
        __DOASSERT( species.end() == find( species.begin(), species.end(), bond.b ), 
                    "Atomic type " << bond.a << " incomplete in input.\n" )

        const std::string bondtype( Physics :: Symbol( bond.b ), Physics :: Symbol( bond.a ) );
        __ASSERT( bonds_.end() != bonds_.find( bondtype ), "Bond already exists.\n" )
        const atomic_specie A( *find( species.begin(), species.end(), bond.a ) );
        const atomic_specie B( *find( species.begin(), species.end(), bond.b ) );
        const types::t_real radius( ( A.radius + B.radius ) * bond.rsigma );
        const types::t_real radius6( radius * radius * radius );
        const types::t_real radius12( radius6 * radius6 );
        LennardJones[ bondtype ].hard_sphere = radius12 * bond.epsilon;
        LennardJones[ bondtype ].van_der_walls = radius6 * bond.epsilon;
      }

      // Hard coded in fortran.
      boost::tuples::get<0>(mesh) = 3;
      boost::tuples::get<1>(mesh) = 3;
      boost::tuples::get<2>(mesh) = 3;

      __TRYEND(, "Could not parse input.\n" )
    }
  } // namespace CLJ
} // namespace LaDa
