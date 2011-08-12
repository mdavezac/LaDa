#ifndef LADA_CRYSTAL_SUPERCELL_H
#define LADA_CRYSTAL_SUPERCELL_H
#include "LaDaConfig.h"

#include <cmath>

#include <Eigen/LU> 

#include <math/misc.h>
#include <math/fuzzy.h>

#include <boost/foreach.hpp>

#include <math/fuzzy.h>
#include <math/smith_normal_form.h>
#include "structure.h"
#include "compare_sites.h"
#include "exceptions.h"
#include "utilities.h"


namespace LaDa
{
  namespace crystal 
  {

    template<class T_TYPE, class T_DERIVED>
      TemplateStructure<T_TYPE> supercell( TemplateStructure<T_TYPE> const &_lattice,
                                           Eigen::DenseBase<T_DERIVED> const &_supercell )
      {
        namespace bt = boost::tuples;
        TemplateStructure<T_TYPE> result; 
        result.cell() = _supercell;

        math::t_SmithTransform transform = math::smith_transform( _lattice.cell(), result.cell());
      
        math::iVector3d &smith = bt::get<1>(transform);
        const math::rMatrix3d factor(bt::get<0>(transform).inverse());
        math::rMatrix3d inv_cell( result.cell().inverse() ); 
        result.reserve(smith(0)*smith(1)*smith(2)*_lattice.size());
        typedef typename TemplateStructure<T_TYPE>::iterator t_iterator;
        typedef typename TemplateStructure<T_TYPE>::const_iterator t_citerator;
        t_citerator const i_site_begin = _lattice.begin();
        t_citerator const i_site_end = _lattice.end();
        
        for( size_t i(0); i < smith(0); ++i )
          for( size_t j(0); j < smith(1); ++j )
            for( size_t k(0); k < smith(2); ++k )
            {
              // in cartesian.
              const math::rVector3d vec( factor * math::rVector3d(i,j,k) );
            
              // adds all lattice sites.
              size_t l(0);
              for( t_citerator i_site(i_site_begin); i_site != i_site_end; ++i_site, ++l)
              {
                Atom<T_TYPE> atom;
                atom.pos    = into_cell(vec+i_site->pos, result.cell(), inv_cell);
                atom.type   = i_site->type;
                atom.freeze = i_site->freeze;
                atom.site   = l;
                result.push_back(atom);
              }
            }
      
        return result;
      }
  } // namespace crystal
} // namespace LaDa
#endif
