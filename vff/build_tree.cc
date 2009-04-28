//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <algorithm>
#include <functional>

#include <opt/smith_normal_form.h>
#include <opt/ndim_iterator.h>
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
      __DOASSERT( structure.lattice->sites.size() != 2,
                  "Lattice should contain 2 different sites.\n" )
      const size_t Nsites( structure.lattice->sites.size() );
      const size_t Natoms( structure.atoms.size() );

      // computes deformation.
      // 16 is the number of first+second neighbors. I think.
      const atat::rMatrix3d deformation( Crystal::retrieve_deformation( structure, 16 ) );

      // now finds smith normal form of ideal lattice.
      atat::iVector3d modulo;
      const atat::rMatrix3d toSmith
      ( 
        to_smith_matrix( deformation, structure.lattice->cell, structure.cell, modulo )
      );
      __ASSERT( modulo(0) * modulo(1) * modulo(2) * 2 != structure.atoms.size(),
                   "Number of atoms in real and deformed matrices do not correspond: "
                <<  modulo(0) << " * " << modulo(1) << " * " <<  modulo(2) 
                << " * " << Nsites << " = " 
                <<  modulo(0) * modulo(1) * modulo(2) * Nsites << " != " 
                <<  structure.atoms.size() << ".\n" )

      // finds ideal first neighbor positions for each lattice sits.
      typedef std::vector< std::vector< atat::rVector3d > > t_FirstNeighbors;
      t_FirstNeighbors first_neighbors( Nsites );
      foreach( const Crystal::Lattice::t_Site &site, structure.lattice->sites )
        first_neighbors[0].push_back( site.pos );
      first_neighbors[1] = first_neighbors[0];
      std::swap( first_neighbors[1][0], first_neighbors[1][1] ); 
      
      for( size_t i(0); i < Nsites; ++i )
      {
        Crystal::find_first_neighbors( first_neighbors[i], structure.lattice->cell, 4 );
        foreach( atat::rVector3d &pos, first_neighbors[i] )
          pos = deformation * pos; // adds in deformation for convenience.
      }
      const types::t_unsigned neighbors_site[2] = { 1, 0 };

      // creates an array indexing each atom.
      types::t_real indices[Nsites][ modulo(0) ][ modulo(1) ][ modulo(2) ];
      for( size_t i(0); i < Nsites; ++i )
        for( size_t j(0); j < modulo(0); ++j )
          for( size_t k(0); k < modulo(1); ++k )
            for( size_t u(0); u < modulo(2); ++u )
              indices[i][j][k][u] = Natoms;
      {
        size_t index(0);
        bool error = false;
        foreach( const Crystal::Structure::t_Atom &atom, structure.atoms )
        {
          atat::iVector3d sindex;
          __ASSERT( atom.site < 0, "site indexing is incorrect.\n" );
          __ASSERT( atom.site > structure.lattice->sites.size(),
                    "site indexing is incorrect.\n" );
          const unsigned site( atom.site );
          smith_index_
          ( 
            toSmith, modulo, 
            atom.pos - structure.lattice->sites[ site ].pos, 
            sindex 
          );
          if( indices[ site ][ sindex(0) ][ sindex(1) ][ sindex(2) ] != Natoms )
          { 
            const size_t i( indices[ site ][ sindex(0) ][ sindex(1) ][ sindex(2) ] );
            std::cerr << "error: " << site << " " << sindex(0) << " " << sindex(1) << " " << sindex(2) << "\n"
                      << i << " " << structure.atoms[i] << "\n" << index << " " << atom << "\n";
            error = true;
          }
          else indices[ site ][ sindex(0) ][ sindex(1) ][ sindex(2) ] = index;
          ++index;
        }
        __ASSERT( error,  "Two sites with same smith index.\n" )
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

      t_Centers :: iterator i_center = centers.begin();
      t_Centers :: iterator i_center_end = centers.end();
      const atat::rMatrix3d inv_cell( !structure.cell );
      i_atom = structure.atoms.begin();
      for(; i_center != i_center_end; ++i_center, ++i_atom )
      {
        const unsigned center_site( i_atom->site );
        const unsigned neighbor_site( neighbors_site[center_site] );
        const atat::rVector3d pos
        ( 
          i_atom->pos - structure.lattice->sites[neighbor_site].pos 
        );
        t_cit i_neigh( first_neighbors[center_site].begin() );
        const t_cit i_neigh_end( first_neighbors[center_site].end() );
        for(; i_neigh != i_neigh_end; ++i_neigh )
        {
          // computes index of nearest neighbor.
          atat::iVector3d sindex;
          smith_index_
          (
            toSmith, modulo, 
            pos + (*i_neigh), 
            sindex
          );
          const types::t_int cindex( indices[neighbor_site][sindex(0)][sindex(1)][sindex(2)] );
          __DOASSERT( cindex == -1, "Index corresponds to no site.\n" )
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
      return true;
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
      for( size_t i(0); i < 3; ++i )
      {
        __DOASSERT
        (
          std::abs( pos(i) - types::t_real( int_pos(i) ) ) > 0.01, 
          "Structure is not ideal.\n"
        )
        _index(i) = int_pos(i) % _modulo(i);
        if( _index(i) < 0 ) _index(i) += _modulo(i);
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

  }
}
