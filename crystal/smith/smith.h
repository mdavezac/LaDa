#ifndef LADA_CRYSTAL_SMITH_H
#define LADA_CRYSTAL_SMITH_H

#include "LaDaConfig.h"

#include <boost/type_traits/is_same.hpp>
#include <boost/static_assert.hpp>


#include <math/smith_normal_form.h>
#include <python/object.h>
#include "../structure/structure.h"
#include "pybase.h"

namespace LaDa 
{
  namespace crystal 
  {
    //! Convenience wrapper around the smuth transform.
    class SmithTransform : public python::Object
    {
        //! \brief Initializer constructor.
        //! \details private so it cannot be constructed without a call throw
        //! smith_transform.
        SmithTransform(PyObject *_in) : python::Object(_in) {}
        //! \brief Initializer constructor.
        //! \details private so it cannot be constructed without a call throw
        //! smith_transform.
        SmithTransform() : python::Object() {}

      public:
        //! Copy constructor.
        SmithTransform(const SmithTransform &_in) : python::Object(_in) {}
        //! Initialization constructor.
        template<class T0, class T1> 
          SmithTransform( Eigen::DenseBase<T0> const &_lattice,
                          Eigen::DenseBase<T1> const &_supercell )
            { init_(_lattice, _supercell); }
        //! Initialization constructor.
        SmithTransform(Structure const &_lattice, Structure const &_supercell)
          { init_(_lattice->cell, _supercell->cell); }

        //! Returns constant reference to transform object.
        math::rMatrix3d const & transform() const 
          { return ((SmithTransformData* const)object_)->transform; }
        //! Returns reference to transform object.
        math::rMatrix3d & transform() 
          { return ((SmithTransformData*)object_)->transform; }
        //! Returns constant reference to quotient object.
        math::iVector3d const & quotient() const 
          { return ((SmithTransformData* const)object_)->quotient; }
        //! Returns reference to quotient object.
        math::iVector3d & quotient() 
          { return ((SmithTransformData*)object_)->quotient; }

#       include "macro.hpp"
        //! Computes smith indices of position \a _pos.
        inline math::iVector3d indices(math::rVector3d const &_pos) const
        {
          LADA_SMITHTRANSFORM_SHARED1(quotient(), transform(), _pos, LADA_PYTHROW,);
          return vector_result;
        }
        //! \brief Computes linear smith index from non-linear smith index.
        inline size_t flat_index(math::iVector3d const &_index, int _site=-1)
        {
          LADA_SMITHTRANSFORM_SHARED0(quotient(), _index, _site);
          return flat_result;
        }
        //! Computes linear smith index of position \a _pos.
        inline size_t flat_index(math::rVector3d const &_pos, int _site=-1)
        {
          LADA_SMITHTRANSFORM_SHARED1(quotient(), transform(), _pos, LADA_PYTHROW,);
          LADA_SMITHTRANSFORM_SHARED0(quotient(), vector_result, _site);
          return flat_result;
        }
        //! Number of unit-cells in the supercell.
        size_t size() const { return LADA_SMITHTRANSFORM_SHARED2(quotient()); }
#       include "macro.hpp"
      private:
        //! creates a smith transform from scratch.
        template<class T0, class T1> 
          void init_(Eigen::DenseBase<T0> const &_lattice, Eigen::DenseBase<T1> const &_supercell);
    };

    template<class T0, class T1> 
      void SmithTransform::init_( Eigen::DenseBase<T0> const &_lattice, 
                                  Eigen::DenseBase<T1> const &_supercell )
      {
        BOOST_STATIC_ASSERT((boost::is_same<typename Eigen::DenseBase<T0>::Scalar,
                                            types::t_real>::value));
        BOOST_STATIC_ASSERT((boost::is_same<typename Eigen::DenseBase<T1>::Scalar,
                                            types::t_real>::value));
        python::Object dummy(object_);
        object_ = smithtransform_type()->tp_alloc(smithtransform_type(), 0);
        if(not object_) return;
        if(not smith_transform_init((SmithTransformData*)object_, _lattice, _supercell)) release();
      }

  } // namespace Crystal

} // namespace LaDa

#endif