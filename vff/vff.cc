#include "LaDaConfig.h"

#include <cstdlib>
#include <algorithm>
#include <functional>
#include <iomanip>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#ifdef _MPI
# include <boost/mpi/collectives/all_reduce.hpp>
#endif

#include <physics/physics.h>
#include <opt/debug.h>
#include <opt/tinyxml.h>
#include <opt/path.h>
#include <crystal/ideal_lattice.h>
#include <mpi/mpi_object.h>


#include "vff.h"
  
namespace LaDa
{
  namespace Vff
  { 

    bool Vff :: Load ( const TiXmlElement &_element )
    {
      namespace bfs = boost::filesystem;

      // some consistency checking
      __ASSERT( not structure.lattice, "Lattice has not been set.\n" )
      if ( structure.lattice->get_nb_sites() != 2 )
      { 
        std::cerr << "Cannot do vff on this lattice.\n" 
                     "Need 2 and only 2 different sites per unit cell.\n";
        return false;
      }
      if (     structure.lattice->get_nb_types(0) + structure.lattice->get_nb_types(1) < 2
           and structure.lattice->get_nb_types(0) + structure.lattice->get_nb_types(1) > 4 )
      { 
        std::cerr << "Cannot do vff on this lattice.\n" 
                     "Need two sites with at most two different atomic types.\n";
        return false;
      }

      const TiXmlElement* parent
        = opt::find_node( _element, "Functional", "type", "vff" );

      if( not parent )
      {
        std::cerr << "Could not find an <Functional type=\"vff\"> tag in input file.\n";
        return false;
      }
      if( not parent->Attribute( "filename" ) ) return load_( *parent );

      const bfs::path path( opt::expand_path( parent->Attribute( "filename" ) ) );
      __DOASSERT( not bfs::exists( path ), path.string() + " does not exist.\n" )
      TiXmlDocument doc;
      opt::read_xmlfile( path, doc );
      __DOASSERT( not doc.FirstChild( "Job" ),
                  "Root tag <Job> does not exist in " + path.string() + ".\n" )
      parent = opt::find_node( *doc.FirstChildElement( "Job" ),
                               "Functional", "type", "vff" );

      if( parent ) return load_( *parent );
      std::cerr << "Could not find an <Functional type=\"vff\"> tag in input file.\n";
      return false;
    }

    bool Vff :: load_( const TiXmlElement &_element )
    {
      std::string str;
      const bool same_species( Crystal::lattice_has_same_species( *structure.lattice ) );


      // reads and initializes bond cutoff
      _element.Attribute( "cutoff", &bond_cutoff );
      if ( bond_cutoff == 0 ) bond_cutoff = types::t_real(1.25); 
      bond_cutoff *= std::sqrt(3.0) / types::t_real(4); // axes bs same as CE
      bond_cutoff *= bond_cutoff; // squared for simplicity
      
      // loads functionals.
      const Crystal :: Lattice :: t_Site &site0 = structure.lattice->sites[0];
      const Crystal :: Lattice :: t_Site &site1 = structure.lattice->sites[1];
      const size_t nb0( site0.type.size() );
      const size_t nb1( site1.type.size() );
      functionals.resize( nb0 + nb1 );
      functionals[0].load( _element, 0u, 0u, site1, structure );
      if( nb0 == 2 ) functionals[1].load( _element, 0u, 1u, site1, structure );
      functionals[nb0].load( _element, 1u, 0u, site0, structure );
      if( nb1 == 2 ) functionals[nb0+1].load( _element, 1u, 1u, site0, structure );

      return true;
    }  // Vff :: Load_


    types::t_real Vff :: energy() const
    {
      types::t_real energy = 0;
      
      LADA_MPI_SPLIT_LOOP( t_Centers :: const_iterator, center, centers, MPI_COMM )
// #     ifdef _MPI
//         for(size_t i(0); i < MPI_COMM.size(); ++i )
//         {
//           if( i != MPI_COMM.rank() ) { MPI_COMM.barrier(); continue; }
//           
//           std::cout << "proc: " << MPI_COMM.rank() << " -- "
//                     << types::t_int( i_center - centers.begin() )
//                     << " -> "
//                     << types::t_int( i_center_end - centers.begin() ) << " "
//                     << "nperproc: " << centernperproc << " "
//                     << "remainder: " << centerremainder << " "
//                     << "b: " << centerbegin_index << " "
//                     << "e: " << centerend_index << " "
//                     << "\n";
//           MPI_COMM.barrier();
//         }
// #     endif
      for (; i_center != i_center_end; ++i_center)
        energy += functionals[i_center->kind()].evaluate( *i_center );
      __MPICODE
      ( 
        energy = boost::mpi::all_reduce( MPI_COMM, energy, std::plus<types::t_real>() ); 
      )

      return energy;
    }

    void Vff :: print_escan_input( const t_Path &_f ) const
    {
      namespace bfs = boost::filesystem;
      std::ostringstream stream;
      types::t_unsigned nb_pseudos=0;

      // prints cell vectors in units of a0 and other
      // whatever nanopes other may be
      for( types::t_unsigned i = 0; i < 3; ++i )
        stream << std::fixed << std::setprecision(7) 
               << std::setw(18) << std::setprecision(8)
               << structure.cell(0,i) * structure.scale / Physics::a0("A")
               << std::setw(18) << std::setprecision(8)
               << structure.cell(1,i) * structure.scale / Physics::a0("A")
               << std::setw(18) << std::setprecision(8)
               << structure.cell(2,i) * structure.scale / Physics::a0("A")
               << std::setw(18) << std::setprecision(8) << structure.cell(0,i) 
               << std::setw(18) << std::setprecision(8) << structure.cell(1,i) 
               << std::setw(18) << std::setprecision(8) << structure.cell(2,i) << "\n";

      // prints atomic position, strain, weight, and atomic position in
      // "other unit"
      t_Centers :: const_iterator i_center = centers.begin();
      t_Centers :: const_iterator i_end = centers.end();
      for(; i_center != i_end; ++i_center )
      {
        // first gets pseudo index
        Crystal::StrAtom stratom;
        __TRYDEBUGCODE( 
          structure.lattice->convert_Atom_to_StrAtom(
             structure.atoms[i_center->get_index()], stratom );, 
             "Error while printing escan input for atom " 
          << i_center->get_index() << ": \n"
          << structure.atoms[i_center->get_index()] << "\n" << structure )
          
        types::t_unsigned index = Physics::Atomic::Z( stratom.type );
        types::t_real msstrain = functionals[i_center->kind()]
                                            .MicroStrain( *i_center, structure );

        // finally goes over bonds and finds number of pseudos and their
        // weights
        AtomicCenter :: const_iterator i_bond = i_center->begin();
        AtomicCenter :: const_iterator i_bond_end = i_center->end();
        typedef std::pair<types::t_unsigned, types::t_unsigned > t_pseudo;
        typedef std::vector< t_pseudo > t_pseudos;
        t_pseudos pseudos;
        for(; i_bond != i_bond_end; ++i_bond )
        { 
          __TRYDEBUGCODE( 
            structure.lattice->convert_Atom_to_StrAtom( 
              structure.atoms[i_bond->get_index()], stratom );, 
               "Error while printing escan input for atoms\n" 
            << structure.atoms[i_bond->get_index()] << "\n" )
          types::t_unsigned Z = Physics::Atomic::Z( stratom.type );
          t_pseudos::iterator i_pseudo = pseudos.begin();
          t_pseudos::iterator i_pseudo_end = pseudos.end();
          for(; i_pseudo != i_pseudo_end; ++i_pseudo )
            if ( i_pseudo->first == Z ) break;

          if ( i_pseudo == i_pseudo_end ) pseudos.push_back( t_pseudo( Z, 1 ) ); 
          else  ++(i_pseudo->second); 
        }

        // now goes over found pseudos and creates output
        t_pseudos::const_iterator i_pseudo = pseudos.begin();
        t_pseudos::const_iterator i_pseudo_end = pseudos.end();
        for( ; i_pseudo != i_pseudo_end; ++i_pseudo )
        {
          math::rVector3d pos = (!structure.cell) * i_center->Origin().pos;
          for( size_t i(0); i < 3; ++i ) pos[i] -= std::floor( pos[i] + 0.000001 );
          ++nb_pseudos;
          stream << std::fixed    << std::setprecision(7)
                 << std::setw(6)  << index << '0'
                    << std::setw(2) << std::setfill('0') 
                    << std::right << i_pseudo->first  
                    << " " << std::setfill(' ') // pseudo index
                 << std::setw(12) << pos[0] << " "  // pseudo position
                 << std::setw(12) << pos[1] << " "  
                 << std::setw(12) << pos[2] << " "  
                 << std::setw(18) << msstrain << " " // microscopic strain
                 << std::setw(6) << std::setprecision(2)
                    << types::t_real( i_pseudo->second ) * 0.25  << " " // weight
                 << std::setw(18) << std::setprecision(7) << pos[0] << " " // pseudo position
                 << std::setw(12) << pos[1] << " "
                 << std::setw(12) << pos[2] << "\n";
        }

      }
      const t_Path directory( _f.parent_path() );
      __TRYBEGIN
        if( not ( directory.empty() or bfs::exists( directory ) ) )
          bfs::create_directory( directory );
        std::ofstream file( _f.string().c_str(), std::ios_base::out|std::ios_base::trunc ); 
        __DOASSERT( file.bad(), "Could not open file " << _f << ".\n" ) 
        // prints number of atoms
        file << nb_pseudos << "\n";
        // print rest of file
        file << stream.str();
        file.flush();
        file.close();
        __ASSERT( not bfs::exists( _f ), _f << " was not created.\n" )
      __TRYEND(, "")
    }

    void Vff :: print_out( std::ostream &stream ) const
    {
      t_AtomicFunctionals :: const_iterator i_func = functionals.begin();
      t_AtomicFunctionals :: const_iterator i_func_end = functionals.end();
      for(; i_func != i_func_end; ++i_func )
        i_func->print_out(stream);
    }

#   ifdef _LADADEBUG
      void Vff :: check_tree() const
      {
        t_Centers :: const_iterator i_center = centers.begin();
        t_Centers :: const_iterator i_center_end = centers.end();
        for(size_t index(0); i_center != i_center_end; ++i_center, ++index )
        {
          __DOASSERT( not i_center->origin, 
                      "Origin of the center is invalid\n"; )
          __DOASSERT( not i_center->structure, 
                      "Invalid pointer to structure\n"; )
          __DOASSERT( i_center->bonds.size() != 4,
                         "Invalid number of bonds: "
                      << i_center->bonds.size() << ", " << index << "\n"; )
          __DOASSERT( i_center->translations.size() != 4,
                         "Invalid number of translations: "
                      << i_center->translations.size() << "\n"; )
        }
      }
#   endif
    boost::tuples::tuple< const types::t_real&, const types::t_real&, const types::t_real&,
                          const types::t_real&, const types::t_real&, const types::t_real& >
      Vff::get_bond( const std::string &_type ) const
      {
        typedef t_AtomicFunctionals :: value_type t_AtomicFunctional;
        namespace bt = boost::tuples;
      
        __DOASSERT( not (structure.lattice and structure.lattice->sites.size() == 2),
                    "Lattice undefined or does not have two sites.\n" )
        namespace bx = boost::xpressive;
        const std::string bond( boost::algorithm::trim_copy( _type ) );
        bx::smatch what;
      
        const bx::sregex regex =(    ( bx::s1 = ( bx::alpha >> !bx::alpha ) )
                                  >> *bx::_s >> "-" >> *bx::_s
                                  >> ( bx::s2 = ( bx::alpha >> !bx::alpha ) ) );
        if( !bx::regex_match( bond, what, regex ) )
          __DOASSERT( true, "Could not find bond " + _type + "\n" )
      
        const bool two_species( structure.lattice->sites[0].type.size() == 2 );
        const std::string A( what.str(1) );
        const std::string B( what.str(2) );
        for( size_t site(0); site < 2; ++site )
        { 
          const int
            i = details::type_index()( structure.lattice->sites[ site ].type,  A ),
            j = details::type_index()( structure.lattice->sites[ site ? 0: 1 ].type,  B );
          if( i == -1 or j == -1 ) continue;
          const size_t Akind
          (
            site == 0 ? size_t(i):( two_species ? size_t(2 + i): size_t(1 + i) )
          );
          const size_t bond_kind
          (
            site == 0 ? size_t(j):( two_species ? size_t(i): 0 )
          );
          __DOASSERT( Akind >= functionals.size(), "Index out-of-range.\n" )
          return functionals[Akind].get_bond( bond_kind );
        }
        __DOASSERT( true, "Could not find angle.\n" )
      }

    boost::tuples::tuple< const types::t_real&, const types::t_real&, 
                          const types::t_real&, const types::t_real&,
                          const types::t_real&, const types::t_real&,
                          const types::t_real& >
      Vff::get_angle( const std::string &_type ) const
      {
        typedef t_AtomicFunctionals :: value_type t_AtomicFunctional;
        namespace bt = boost::tuples;
      
        __DOASSERT( not (structure.lattice and structure.lattice->sites.size() == 2),
                    "Lattice undefined or does not have two sites.\n" )
        namespace bx = boost::xpressive;
        const std::string bond( boost::algorithm::trim_copy( _type ) );
        bx::smatch what;
      
        const bx::sregex regex =(    ( bx::s1 = ( bx::alpha >> !bx::alpha ) )
                                  >> *bx::_s >> "-" >> *bx::_s
                                  >> ( bx::s2 = ( bx::alpha >> !bx::alpha ) )
                                  >> *bx::_s >> "-" >> *bx::_s
                                  >> ( bx::s3 = ( bx::alpha >> !bx::alpha ) ) );
        if( !bx::regex_match( bond, what, regex ) )
          __DOASSERT( true, "Could not find angle " + _type + "\n" )
      
        const bool two_species( structure.lattice->sites[0].type.size() == 2 );
        const std::string A( what.str(1) );
        const std::string B( what.str(2) );
        const std::string C( what.str(3) );
        for( size_t site(1); site < 2; ++site )
        { 
          const int
            i = details::type_index()( structure.lattice->sites[ site ? 0: 1 ].type,  A ),
            j = details::type_index()( structure.lattice->sites[ site ].type,  B ),
            k = details::type_index()( structure.lattice->sites[ site ? 0: 1 ].type,  C );
          if( i == -1 or j == -1 or k == -1 ) continue;
          const size_t Bkind
          (
            site == 0 ? size_t(j):( two_species ? size_t(2 + j): size_t(1 + j) )
          );
          const size_t angle_kind
          (
              size_t( site == 0 ? size_t(i):( two_species ? size_t(i): 0 ) )
            + size_t( site == 0 ? size_t(k):( two_species ? size_t(k): 0 ) )
          );
          __DOASSERT( Bkind >= functionals.size(), "Index out-of-range.\n" )
          return functionals[Bkind].get_angle( angle_kind );
        }
        __DOASSERT( true, "Could not find angle.\n" )
      }
  } // namespace vff
} // namespace LaDa
