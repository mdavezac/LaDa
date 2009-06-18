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
#include "smith.h"
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
          namespace bt = boost::tuples;
          __ASSERT( _structure.lattice == NULL, "Lattice not set.\n" )
          t_SmithTransform transform
            = get_smith_transform( _structure.lattice->cell, _structure.cell);
        
          atat::iVector3d &smith = bt::get<1>(transform);
          const atat::rMatrix3d factor
          ( 
             (!_structure.cell) * (!bt::get<0>(transform))
          );
          typename TStructure<T_TYPE>::t_Atom atom;
          for( size_t i(0); i < smith(0); ++i )
            for( size_t j(0); j < smith(1); ++j )
              for( size_t k(0); k < smith(2); ++k )
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
                const atat::rVector3d vec( _structure.cell * vec2 );
              
                // adds all lattice sites.
                typedef Crystal::Lattice::t_Site t_Site;
                size_t i(0);
                foreach( const t_Site &site, _structure.lattice->sites ) 
                {
                  atom.pos = vec + site.pos;
                  atom.site = i;
                  atom.freeze = site.freeze;
                  _structure.atoms.push_back(atom);
                  ++i;
                }
              }
        
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
      bool const result( details::fill_structure( _str ) );
      if( not result ) return false;
      foreach( Crystal::TStructure<std::string>::t_Atom & atom, _str.atoms )
        atom.type = _str.lattice->sites[ atom.site ].type[0];
      return result;
    }

  } // namespace Crystal

} // namespace LaDa
