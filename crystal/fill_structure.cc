#include "LaDaConfig.h"

#include <functional>

#include <math/smith_normal_form.h>
#include <opt/debug.h>

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
          LADA_NASSERT( _structure.lattice == NULL, "Lattice not set.\n" )
          t_SmithTransform transform
            = get_smith_transform( _structure.lattice->cell, _structure.cell);
        
          math::iVector3d &smith = bt::get<1>(transform);
          const math::rMatrix3d factor
          ( 
             (!bt::get<0>(transform))
          );
          math::rMatrix3d inv_cell( _structure.cell.inverse() ); 
          _structure.atoms.resize(smith(0)*smith(1)*smith(2)*_structure.lattice->sites.size());
          typedef typename Crystal::TStructure<T_TYPE>::t_Atoms::iterator t_iterator;
          t_iterator i_atom = _structure.atoms.begin();
#         ifdef LADA_DEBUG
             t_iterator const i_atom_end = _structure.atoms.end();
#         endif
          Crystal::Lattice::t_Sites::iterator const i_site_begin = _structure.lattice->sites.begin();
          Crystal::Lattice::t_Sites::iterator const i_site_end = _structure.lattice->sites.end();
          
          for( size_t i(0); i < smith(0); ++i )
            for( size_t j(0); j < smith(1); ++j )
              for( size_t k(0); k < smith(2); ++k )
              {
                // in cartesian.
                const math::rVector3d vec( factor * math::rVector3d(i,j,k) );
              
                // adds all lattice sites.
                typedef Crystal::Lattice::t_Site t_Site;
                size_t l(0);
                for( Crystal::Lattice::t_Sites::iterator i_site(i_site_begin);
                     i_site != i_site_end; ++i_site, ++i_atom )
                {
                  LADA_ASSERT(i_atom != i_atom_end, "out of range.")
                  // in supercell fractional.
                  math::rVector3d frac( inv_cell * ( vec + i_site->pos ) );
                  // in supercell fractional and in supercell parallelogram
                  const math::rVector3d inside
                  (
                    frac(0) - std::floor( frac(0) + 0.0000001 ),
                    frac(1) - std::floor( frac(1) + 0.0000001 ),
                    frac(2) - std::floor( frac(2) + 0.0000001 )
                  );
                  // back to cartesian.
                  i_atom->pos = _structure.cell * inside;
                  i_atom->site = l;
                  i_atom->freeze = i_site->freeze;
                  ++l;
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
