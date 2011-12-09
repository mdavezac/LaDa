#ifndef LADA_CRYSTAL_SUPERCELL_H
#define LADA_CRYSTAL_SUPERCELL_H
#include "LaDaConfig.h"

#include <cmath>

#include <Eigen/LU> 

#include <math/misc.h>
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
      Structure<T_TYPE> supercell( Structure<T_TYPE> const &_lattice,
                                           Eigen::DenseBase<T_DERIVED> const &_supercell )
      {
        namespace bt = boost::tuples;
        Structure<T_TYPE> result; 
        result->cell = _supercell;
        result->scale = _lattice->scale;
        if(_lattice->name.size() != 0) result->name = "supercell of " + _lattice->name;
        result->energy = _lattice->energy * types::t_real(result.size())/types::t_real(_lattice.size());

        math::t_SmithTransform transform = math::smith_transform( _lattice.cell(), result.cell());
      
        math::iVector3d &smith = bt::get<1>(transform);
        const math::rMatrix3d factor(bt::get<0>(transform).inverse());
        math::rMatrix3d inv_cell( result.cell().inverse() ); 
        result.reserve(smith(0)*smith(1)*smith(2)*_lattice.size());
        typedef typename Structure<T_TYPE>::iterator t_iterator;
        typedef typename Structure<T_TYPE>::const_iterator t_citerator;
        t_citerator const i_site_begin = _lattice.begin();
        t_citerator const i_site_end = _lattice.end();
        
        for( math::iVector3d::Scalar i(0); i < smith(0); ++i )
          for( math::iVector3d::Scalar j(0); j < smith(1); ++j )
            for( math::iVector3d::Scalar k(0); k < smith(2); ++k )
            {
              // in cartesian.
              const math::rVector3d vec( factor * math::rVector3d(i,j,k) );
            
              // adds all lattice sites.
              types::t_unsigned l(0);
              for( t_citerator i_site(i_site_begin); i_site != i_site_end; ++i_site, ++l)
              {
                Atom<T_TYPE> atom( into_cell(vec+i_site->pos(), result->cell, inv_cell),
                                   i_site->type(), l, i_site->freeze() );
                result.push_back(atom);
              }
            }
      
        return result;
      }
  } // namespace crystal
} // namespace LaDa
#endif
