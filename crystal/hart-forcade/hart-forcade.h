#ifndef LADA_CRYSTAL_HF_H
#define LADA_CRYSTAL_HF_H

#include "LaDaConfig.h"

#include <boost/type_traits/is_same.hpp>
#include <boost/static_assert.hpp>


#include <math/misc.h>
#include <math/smith_normal_form.h>
#include <python/object.h>
#include "../structure/structure.h"
#include "pybase.h"

namespace LaDa 
{
  namespace crystal 
  {
    //! Convenience wrapper around the smuth transform.
    class HFTransform : public python::Object
    {
        //! \brief Initializer constructor.
        //! \details private so it cannot be constructed without a call throw
        //! hf_transform.
        HFTransform(PyObject *_in) : python::Object(_in) {}
        //! \brief Initializer constructor.
        //! \details private so it cannot be constructed without a call throw
        //! hf_transform.
        HFTransform() : python::Object() {}

      public:
        //! Copy constructor.
        HFTransform(const HFTransform &_in) : python::Object(_in) {}
        //! Initialization constructor.
        template<class T0, class T1> 
          HFTransform( Eigen::DenseBase<T0> const &_lattice,
                          Eigen::DenseBase<T1> const &_supercell )
            { init_(_lattice, _supercell); }
        //! Initialization constructor.
        HFTransform(Structure const &_lattice, Structure const &_supercell)
          { init_(_lattice->cell, _supercell->cell); }

        //! Returns constant reference to transform object.
        math::rMatrix3d const & transform() const 
          { return ((PyHFTObject* const)object_)->transform; }
        //! Returns reference to transform object.
        math::rMatrix3d & transform() 
          { return ((PyHFTObject*)object_)->transform; }
        //! Returns constant reference to quotient object.
        math::iVector3d const & quotient() const 
          { return ((PyHFTObject* const)object_)->quotient; }
        //! Returns reference to quotient object.
        math::iVector3d & quotient() 
          { return ((PyHFTObject*)object_)->quotient; }

#       include "macro.hpp"
        //! Computes hf indices of position \a _pos.
        inline math::iVector3d indices(math::rVector3d const &_pos) const
        {
          LADA_HFTRANSFORM_SHARED1(quotient(), transform(), _pos, LADA_PYTHROW,);
          return vector_result;
        }
        //! \brief Computes linear hf index from non-linear hf index.
        inline size_t flat_index(math::iVector3d const &_index, int _site=-1)
        {
          LADA_HFTRANSFORM_SHARED0(quotient(), _index, _site);
          return flat_result;
        }
        //! Computes linear hf index of position \a _pos.
        inline size_t flat_index(math::rVector3d const &_pos, int _site=-1)
        {
          LADA_HFTRANSFORM_SHARED1(quotient(), transform(), _pos, LADA_PYTHROW,);
          LADA_HFTRANSFORM_SHARED0(quotient(), vector_result, _site);
          return flat_result;
        }
        //! Number of unit-cells in the supercell.
        size_t size() const { return LADA_HFTRANSFORM_SHARED2(quotient()); }
#       include "macro.hpp"
      private:
        //! creates a hf transform from scratch.
        template<class T0, class T1> 
          void init_(Eigen::DenseBase<T0> const &_lattice, Eigen::DenseBase<T1> const &_supercell);
    };

    template<class T0, class T1> 
      void HFTransform::init_( Eigen::DenseBase<T0> const &_lattice, 
                                  Eigen::DenseBase<T1> const &_supercell )
      {
        BOOST_STATIC_ASSERT((boost::is_same<typename Eigen::DenseBase<T0>::Scalar,
                                            types::t_real>::value));
        BOOST_STATIC_ASSERT((boost::is_same<typename Eigen::DenseBase<T1>::Scalar,
                                            types::t_real>::value));
        python::Object dummy(object_);
        object_ = hftransform_type()->tp_alloc(hftransform_type(), 0);
        if(not object_) return;
        if(not _init_hft((PyHFTObject*)object_, _lattice, _supercell)) release();
      }

  } // namespace Crystal

} // namespace LaDa

#endif
