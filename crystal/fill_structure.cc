//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <functional>

#include <opt/ndim_iterator.h>
#include <opt/smith_normal_form.h>
#include <opt/debug.h>
#include <atat/misc.h>

#include "lattice.h"
#include "structure.h"
#include "fill_structure.h"



namespace LaDa
{

  namespace Crystal 
  {
    namespace details
    {
      template< class T_TYPE >
        bool fill_structure( Crystal::TStructure<T_TYPE> &_structure )
        {
          std::cout << " 1.1 \n";
          __ASSERT( _structure.lattice == NULL, "Lattice not set.\n" )
          std::cout << _structure.cell << "\n";
          std::cout << Crystal::Structure::lattice << "\n";
          std::cout << Crystal::TStructure<T_TYPE>::lattice << "\n";
          std::cout << " 1.11 \n";
          const atat::rMatrix3d inv_lat_cell( ( !_structure.lattice->cell ) * _structure.cell );
          std::cout << " 1.12 \n";
          atat::iMatrix3d int_cell;
          for( size_t i(0); i < 3; ++i )
            for( size_t j(0); j < 3; ++j )
            {
          std::cout << " 1.13 " << i << " " << j << "\n";
              const types::t_real d( std::floor( inv_lat_cell(i,j) + 0.5  ) );
              int_cell(i,j) = types::t_int( d );
              __DOASSERT
              ( 
                std::abs( types::t_real( int_cell(i,j) ) - inv_lat_cell(i,j) ) > 0.01,
                   "Input structure is not supercell of the lattice: " 
                << int_cell(i,j) << " != " << inv_lat_cell(i,j) << "\n"
              )
            }
          std::cout << " 1.2 \n";
          atat::iMatrix3d left, right, smith;
          opt::smith_normal_form( smith, left, int_cell, right );
          std::cout << " 1.3 \n";
        
          atat::rMatrix3d rleft;
          for( size_t i(0); i < 3; ++i )
            for( size_t j(0); j < 3; ++j )
              rleft(i,j) = types::t_real( left(i,j) );
          const atat::rMatrix3d factor
          ( 
            ( !_structure.cell ) * _structure.lattice->cell * (!rleft) 
          );
          std::cout << " 1.4 \n";
          typename TStructure<T_TYPE>::t_Atom atom;
          for( size_t i(0); i < smith(0,0); ++i )
            for( size_t j(0); j < smith(1,1); ++j )
              for( size_t k(0); k < smith(2,2); ++k )
              {
                // in supercell fractional
                const atat::rVector3d vec1( factor * atat::rVector3d(i,j,k) );
                // in supercell fractional and in supercell parallelogram
                const atat::rVector3d vec2
                (
                  vec1(0) - std::floor( vec1(0) ),
                  vec1(1) - std::floor( vec1(1) ),
                  vec1(2) - std::floor( vec1(2) )
                );
                // in cartesian
                const atat::rVector3d vec( _structure.cell * vec1 );
              
                // adds all lattice sites.
                typedef Crystal::Lattice::t_Site t_Site;
                foreach( const t_Site &site, _structure.lattice->sites ) 
                {
                  atom.pos = vec + site.pos;
                  atom.site = site.site;
                  atom.freeze = site.freeze;
                  _structure.atoms.push_back(atom);
                }
              }
          std::cout << " 1.5 \n";
        
          return true;
        }
    } 


    bool fill_structure( Crystal::Structure &_str )
    {
      bool const result( details::fill_structure( _str ) );
      if( not result ) return false;
      foreach( Crystal::Structure::t_Atom & atom, _str.atoms ) atom.type = -1e0;
      return result;
    }
   
    bool fill_structure( Crystal::TStructure<std::string> &_str )
    {
      std::cout << " 1 \n";
      std::cout << _str.lattice << "\n";
      bool const result( details::fill_structure( _str ) );
      std::cout << " 2 \n";
      if( not result ) return false;
      std::cout << " 3 \n";
      foreach( Crystal::TStructure<std::string>::t_Atom & atom, _str.atoms )
        atom.type = _str.lattice->sites[ atom.site ].type[0];
      std::cout << " 4 \n";
      return result;
    }

  } // namespace Crystal

} // namespace LaDa
