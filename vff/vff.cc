//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <cstdlib>
#include <algorithm>
#include <functional>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <physics/physics.h>
#include <opt/ndim_iterator.h>
#include <opt/atat.h>
#include <opt/debug.h>
#include <opt/atat.h>
#include <opt/tinyxml.h>
#include <opt/smith_normal_form.h>
#include <crystal/ideal_lattice.h>

#include "vff.h"
  
namespace LaDa
{
  namespace Vff
  { 
    bool Vff :: build_tree()
    {
      __TRYBEGIN
      __DOASSERT( structure.lattice == NULL, "Lattice not set.\n" )

      // computes deformation.
      const atat::rMatrix3d deformation( Crystal::retrieve_deformation( structure, 16 ) );

      // now finds smith normal form of ideal lattice.
      atat::iVector3d modulo;
      const atat::rMatrix3d toSmith
      ( 
        to_smith_matrix( deformation, structure.lattice->cell, structure.cell, modulo )
      );

      // finds ideal first neighbor positions.
      std::vector< atat::rVector3d > first_neighbors;
      foreach( const Crystal::Lattice::t_Site &site, structure.lattice->sites )
        first_neighbors.push_back( site.pos );
      Crystal::find_first_neighbors( first_neighbors, structure.lattice->cell, 4 );
      foreach( atat::rVector3d &pos, first_neighbors )
      {
        pos -= structure.lattice->sites[0].pos;
        pos = deformation * pos; // adds in deformation for convenience.
      }

      // creates an array indexing each atom.
      types::t_real indices[structure.lattice->sites.size()]
                           [ modulo(0) ][ modulo(1) ][ modulo(2) ];
      {
        size_t index(0);
        foreach( const Crystal::Structure::t_Atom &atom, structure.atoms )
        {
          atat::iVector3d sindex;
          __ASSERT( atom.site < 0, "site indexing is incorrect.\n" );
          __ASSERT( atom.site > structure.lattice->sites.size(),
                    "site indexing is incorrect.\n" );
          smith_index_
          ( 
            toSmith, modulo, 
            atom.pos - structure.lattice->sites[ atom.site ].pos, 
            sindex 
          );
          indices[ unsigned(atom.site) ][ sindex(0) ][ sindex(1) ][ sindex(2) ] = index;
          ++index;
        }
      }
      
      // constructs list of centers.
      centers.clear();
      centers.reserve( structure.atoms.size() );
      Crystal::Structure::t_Atoms::iterator i_atom = structure.atoms.begin();
      Crystal::Structure::t_Atoms::iterator i_atom_end = structure.atoms.end();
      for(types::t_unsigned index=0; i_atom != i_atom_end; ++i_atom, ++index )
        centers.push_back( AtomicCenter( structure, *i_atom, index ) );

      // finally builds tree.
      typedef std::vector< atat::rVector3d > :: const_iterator t_cit;
      const t_cit i_neigh_begin( first_neighbors.begin() );
      const t_cit i_neigh_end( first_neighbors.end() );

      t_Centers :: iterator i_center = centers.begin();
      t_Centers :: iterator i_center_end = centers.end();
      const atat::rMatrix3d inv_cell( !structure.cell );
      i_atom = structure.atoms.begin();
      for(; i_center != i_center_end; ++i_center, ++i_atom )
      {
        const unsigned site( i_atom->site == 0 ? 0: 1 );
        const atat::rVector3d pos
        ( 
          i_atom->pos - structure.lattice->sites[site].pos 
        );
        for( t_cit i_neigh( i_neigh_begin ); i_neigh != i_neigh_end; ++i_neigh )
        {
          // computes index of nearest neighbor.
          atat::iVector3d sindex;
          smith_index_
          (
            toSmith, modulo, 
            pos + (*i_neigh), 
            sindex
          );
          const types::t_int cindex( indices[site][ sindex(0) ][ sindex(1) ][ sindex(2) ] );
          // now creates branch in tree.
          t_Centers :: iterator i_bond( centers.begin() + cindex );
          i_center->bonds.push_back( t_Center::__make__iterator__( i_bond ) );
          const atat::rVector3d dfrac
          ( 
              inv_cell 
            * ( 
                  (const atat::rVector3d) *i_center 
                - (const atat::rVector3d) *i_bond
              )
           ); 
          const atat::rVector3d frac
          (
            dfrac(0) - rint( dfrac(0) ),
            dfrac(1) - rint( dfrac(1) ),
            dfrac(2) - rint( dfrac(2) )
          );
          i_center->translations.push_back( frac );
          i_center->do_translates.push_back
          ( 
            atat::norm2(dfrac) > atat::zero_tolerance 
          );
        }
      }
      __ENDGROUP__
      catch( ... )
      {
        std::cerr << "Could not build tree.\n";
        return false;
      }
    }

    void Vff :: smith_index_( const atat::rMatrix3d &_toSmith,
                              const atat::iVector3d &_modulo,
                              const atat::rVector3d &_pos,
                              atat::iVector3d &_index )
    {
      const atat::rVector3d pos( _toSmith * _pos );
      const atat::iVector3d int_pos
      (
        types::t_int( rint( pos(0) ) ),
        types::t_int( rint( pos(1) ) ),
        types::t_int( rint( pos(2) ) )
      );
      std::cout << _pos << ", " << pos << ", " << int_pos << "\n";
      for( size_t i(0); i < 3; ++i )
      {
        __DOASSERT
        (
          std::abs( pos(i) - types::t_real( int_pos(i) ) ) > 0.01, 
          "Structure is not ideal.\n"
        )
        _index(i) = int_pos(i) % _modulo(i);
      }
    }

    atat::rMatrix3d Vff :: to_smith_matrix( const atat::rMatrix3d &_deformation,
                                            const atat::rMatrix3d &_lat_cell,
                                            const atat::rMatrix3d &_str_cell,
                                            atat::iVector3d &_modulo )
    {
      atat::rMatrix3d result;
      atat::iMatrix3d left, right, smith;
      const atat::rMatrix3d inv_lat( !_lat_cell );
      const atat::rMatrix3d inv_lat_cell( inv_lat * _deformation * _str_cell );
      std::cout << _lat_cell << "\n\n" << inv_lat << "\n\n";
      std::cout << _str_cell << "\n\n" << inv_lat_cell << "\n";
      atat::iMatrix3d int_cell;
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
        {
          int_cell(i,j) = types::t_int( rint( inv_lat_cell(i,j) ) ); 
          __DOASSERT
          ( 
            std::abs( types::t_real( int_cell(i,j) ) - inv_lat_cell(i,j) ) > 0.01,
               "Input structure is not supercell of the lattice: " 
            << int_cell(i,j) << " != " << inv_lat_cell(i,j) << "\n"
          )
        }
      opt::smith_normal_form( smith, left, int_cell, right );
      for( size_t i(0); i < 3; ++i )
      {
        for( size_t j(0); j < 3; ++j )
          result(i,j) = types::t_real( left(i,j) );
        _modulo(i) = smith(i,i);
      }
      return result * ( !_lat_cell ) * _deformation;
    }

    bool Vff :: initialize_centers()
    {
      centers.clear();
      
      // Creates a list of centers
      t_Atoms :: iterator i_atom = structure.atoms.begin();
      t_Atoms :: iterator i_atom_end = structure.atoms.end();

      for(types::t_unsigned index=0; i_atom != i_atom_end; ++i_atom, ++index )
        centers.push_back( AtomicCenter( structure, *i_atom, index ) );

      // Creates a list of closest neighbors
      std::vector< atat::rVector3d > neighbors;
      typedef Crystal::Lattice::t_Sites :: iterator t_it;
      t_it i_site_begin = structure.lattice->sites.begin();
      t_it i_site, i_site2;
      t_it i_site_end = structure.lattice->sites.end();
      
      for(i_site = i_site_begin; i_site != i_site_end; ++i_site )
      {
        for(i_site2 = i_site_begin; i_site2 != i_site_end; ++i_site2 )
        {
          opt::Ndim_Iterator< types::t_int, std::less_equal<types::t_int> > period;
          period.add(-1,1);
          period.add(-1,1);
          period.add(-1,1);
          do // goes over periodic image
          {
            // constructs perdiodic image of atom *i_bond
            atat::rVector3d frac_image, image;
            frac_image[0] =  (types::t_real) period.access(0);
            frac_image[1] =  (types::t_real) period.access(1);
            frac_image[2] =  (types::t_real) period.access(2);
            image = i_site2->pos + structure.lattice->cell * frac_image;
            if( atat::norm2( image - i_site->pos ) > types::tolerance )
              neighbors.push_back( image - i_site->pos );
          } while ( ++period ); 
        }
      }

      // Sorts the neighbors according to distance from origin
      std::sort( neighbors.begin(), neighbors.end(), atat::norm_compare() );
      // And reduces to first neighbors only
      neighbors.resize(4*structure.lattice->sites.size());

      t_Centers :: iterator i_begin = centers.begin();
      t_Centers :: iterator i_end = centers.end();
      t_Centers :: iterator i_center, i_bond;
      atat::rVector3d image;
      atat::rVector3d frac_image;
      atat::rVector3d cut;
      cut = neighbors.front();
      types::t_real cutoff = types::t_real(0.25) * atat::norm2( neighbors.front() );
      for( i_center = i_begin; i_center != i_end; ++i_center )
      {
        for( i_bond = i_begin; i_bond != i_end; ++i_bond)
          if( i_bond != i_center )
          {
            std::vector<atat::rVector3d> :: const_iterator i_neigh = neighbors.begin();
            std::vector<atat::rVector3d> :: const_iterator i_neigh_end = neighbors.end();
            for(; i_neigh != i_neigh_end; ++i_neigh )
            {
              image = i_center->origin->pos - *i_neigh - i_bond->origin->pos;
              frac_image = (!structure.cell) * image;
              cut[0] = frac_image[0] - rint( frac_image[0] );
              cut[1] = frac_image[1] - rint( frac_image[1] );
              cut[2] = frac_image[2] - rint( frac_image[2] );
              cut = structure.cell * cut;
              if( atat::norm2( cut ) < cutoff )
              {
                i_center->bonds.push_back( t_Center ::__make__iterator__(i_bond) );
                frac_image[0] = rint( frac_image[0] );
                frac_image[1] = rint( frac_image[1] );
                frac_image[2] = rint( frac_image[2] );
                i_center->translations.push_back( frac_image );
                i_center->do_translates.push_back
                ( 
                  atat::norm2(frac_image) > atat::zero_tolerance 
                );
              }
            }
          }
      }

      __DODEBUGCODE( check_tree(); )
      return true;
    } // Vff :: construct_bond_list

    bool Vff :: construct_centers()
    {
      centers.clear();
      
      t_Atoms :: iterator i_atom = structure.atoms.begin();
      t_Atoms :: iterator i_atom_end = structure.atoms.end();

      for(types::t_unsigned index=0; i_atom != i_atom_end; ++i_atom, ++index )
        centers.push_back( AtomicCenter( structure, *i_atom, index ) );

      t_Centers :: iterator i_begin = centers.begin();
      t_Centers :: iterator i_end = centers.end();
      t_Centers :: iterator i_center, i_bond;

      for( i_center = i_begin; i_center != i_end; ++i_center )
        for( i_bond = i_begin; i_bond != i_end; ++i_bond )
          if ( i_bond != i_center )
            i_center->add_bond ( t_Center::__make__iterator__(i_bond), bond_cutoff);

      
      // consistency check
      for( i_center = i_begin; i_center != i_end; ++i_center )
        if ( i_center->size() != 4 )
        {
          std::cerr << " Atomic center at " << (atat::rVector3d) i_center->Origin()
                    << " has " << i_center->size() 
                    << " bonds!!" << std::endl;
          return false;
        } 

      __DODEBUGCODE( check_tree(); )
      return true;
    } // Vff :: construct_bond_list

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
                     "Need two sites with at most two different atomic types\n";
        return false;
      }

      const TiXmlElement* parent
        = opt::find_node( _element, "Functional", "type", "vff" );

      if( not parent )
      {
        std::cerr << "Could not find an <Functional type=\"vff\"> tag in input file" 
                  << std::endl;
        return false;
      }
      if( not parent->Attribute( "filename" ) ) return load_( *parent );

      const bfs::path path( Print::reformat_home( parent->Attribute( "filename" ) ) );
      __DOASSERT( not bfs::exists( path ), path.string() + " does not exist.\n" )
      TiXmlDocument doc;
      opt::read_xmlfile( path, doc );
      __DOASSERT( not doc.FirstChild( "Job" ),
                  "Root tag <Job> does not exist in " + path.string() + ".\n" )
      parent = opt::find_node( *doc.FirstChildElement( "Job" ), "Functional", "type", "vff" );

      if( parent ) return load_( *parent );
      std::cerr << "Could not find an <Functional type=\"vff\"> tag in input file" 
                << std::endl;
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
      
      t_Centers :: const_iterator i_center = centers.begin();
      t_Centers :: const_iterator i_end = centers.end();
      for (; i_center != i_end; ++i_center)
        energy += functionals[i_center->kind()].evaluate( *i_center );

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
          atat::rVector3d pos = (!structure.cell) * i_center->Origin().pos;
          ++nb_pseudos;
          stream << std::fixed    << std::setprecision(7)
                 << std::setw(6)  << index << '0' << i_pseudo->first  // pseudo index
                 << std::setw(12) << pos[0] // pseudo position
                 << std::setw(12) << pos[1] 
                 << std::setw(12) << pos[2] 
                 << std::setw(18) << msstrain << " " // microscopic strain
                 << std::setw(6) << std::setprecision(2)
                                 << types::t_real( i_pseudo->second ) * 0.25  // weight
                 << std::setw(18) << std::setprecision(7) << pos[0] // pseudo position
                 << std::setw(12) << pos[1] 
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
        for(; i_center != i_center_end; ++i_center )
        {
          __DOASSERT( not i_center->origin, 
                      "Origin of the center is invalid\n"; )
          __DOASSERT( not i_center->structure, 
                      "Invalid pointer to structure\n"; )
          __DOASSERT( i_center->bonds.size() != 4,
                         "Invalid number of bonds: "
                      << i_center->bonds.size() << "\n"; )
          __DOASSERT( i_center->translations.size() != 4,
                         "Invalid number of translations: "
                      << i_center->translations.size() << "\n"; )
        }
      }
#   endif
  } // namespace vff
} // namespace LaDa
