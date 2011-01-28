#include "LaDaConfig.h"

#include <boost/python/def.hpp>
#include <boost/python/object.hpp>
#ifdef LADA_MPI
# include <boost/mpi/python.hpp>
#endif

#include <opt/types.h>
#include <opt/debug.h>
#include <math/serialize.h>
#include <python/numpy_types.h>

#include "../structure.h"
#include "misc.hpp"


namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    namespace mn = math::numpy;
    //! Calls functor _op on vector, as given by last axis.
    template<class T_OP> 
      bp::object unloop_array_call(bp::object const &_array, T_OP &_op)
      {
        if(_array.ptr() == Py_None or _array.ptr() == NULL) 
        { 
          PyErr_SetString(PyExc_ValueError, "Got None where a numpy array is expected.");
          bp::throw_error_already_set();
          return bp::object();
        }
        PyArray_Check(_array.ptr());
        mn::check_is_array(_array);
        PyArrayObject *ptr_array = mn::get_pyarray_pointer(_array);
        if(PyArray_DIM(_array.ptr(), ptr_array->nd-1) != 3)
        {
          PyErr_SetString(PyExc_ValueError, "Last dimension of array should be 3.");
          bp::throw_error_already_set();
          return bp::object();
        }
  
        int last = std::max(ptr_array->nd - 1, 0);
        PyObject *i_in  = PyArray_IterAllButAxis(_array.ptr(), &last);
        int const in_stride = mn::get_strides_pointer(_array)[last];
        int const in_stride2 = in_stride << 1;
        int const type = PyArray_TYPE(_array.ptr());
        _op.init(_array);
        while( PyArray_ITER_NOTDONE(i_in) )
        {
          math::rVector3d const orig( 
             mn::to_type<types::t_real>(i_in, type),
             mn::to_type<types::t_real>(i_in, type, in_stride),
             mn::to_type<types::t_real>(i_in, type, in_stride2) );
          _op(orig);
          PyArray_ITER_NEXT(i_in); 
          ++_op;
        }
        Py_DECREF(i_in);
        return _op.result;
      }

    struct ToCell
    {
      math::rVector3d center;
      math::rMatrix3d cell, invcell;
      boost::python::object iterator;
      boost::python::object result;
      int last;
      int stride;
      int stride2;
      ToCell   (math::rMatrix3d const &_cell, math::rVector3d const &_center)
             : center(_center), cell(_cell) { invcell = cell.inverse(); }
      void operator()(math::rVector3d const &_orig) const
      {
        math::rVector3d const fractional(invcell * (_orig-center));
        math::rVector3d const centered(
            fractional[0] - std::floor(fractional[0] + 0.5),
            fractional[1] - std::floor(fractional[1] + 0.5),
            fractional[2] - std::floor(fractional[2] + 0.5) );
        math::rVector3d const O1( 
            centered[0] > 0e0 ? centered[0] + 1e0: centered[0], 
            centered[1] > 0e0 ? centered[1] + 1e0: centered[1], 
            centered[2] > 0e0 ? centered[2] + 1e0: centered[2] );
        math::rVector3d vec(cell * O1);
        char const * const ptr_out = (char const * const)(PyArray_ITER_DATA(iterator.ptr()));
        *(types::t_real*)(ptr_out)           = vec[0];
        *(types::t_real*)(ptr_out + stride)  = vec[1];
        *(types::t_real*)(ptr_out + stride2) = vec[2];
      }
      void init(bp::object _array)
      {
        result = mn::from_template(_array, mn::type<types::t_real>::value);
        PyArrayObject const * const ptr_array = mn::get_pyarray_pointer(_array);
        last = std::max(ptr_array->nd - 1, 0);
        stride = mn::get_strides_pointer(result)[last];
        stride2 = stride << 1;
        iterator = bp::object( bp::handle<>(PyArray_IterAllButAxis(result.ptr(), &last)));
      }
      void operator++() { PyArray_ITER_NEXT(iterator.ptr()); }
    }; 

    struct ToOrigin : public ToCell
    {
      ToOrigin   (math::rMatrix3d const &_cell, math::rVector3d const &_center)
               : ToCell(_cell, _center) {}
      void operator()(math::rVector3d const &_orig) const
      {
        math::rVector3d const fractional(invcell * (_orig-center));
        math::rVector3d const centered_frac(
            fractional[0] - std::floor(fractional[0] + 0.5),
            fractional[1] - std::floor(fractional[1] + 0.5),
            fractional[2] - std::floor(fractional[2] + 0.5) );
        math::rVector3d const vec(cell * centered_frac);
        char const * const ptr_out = (char const * const)(PyArray_ITER_DATA(iterator.ptr()));
        *(types::t_real*)(ptr_out)           = vec[0];
        *(types::t_real*)(ptr_out + stride)  = vec[1];
        *(types::t_real*)(ptr_out + stride2) = vec[2];
      }
    }; 

    struct ToVoronoi : public ToCell
    {
      ToVoronoi   (math::rMatrix3d const &_cell, math::rVector3d const &_center)
                : ToCell(_cell, _center) {}
      void operator()(math::rVector3d const &_orig) const
      {
        math::rVector3d const fractional(invcell * (_orig-center));
        math::rVector3d const centered_frac(
            fractional[0] - std::floor(fractional[0] + 0.5),
            fractional[1] - std::floor(fractional[1] + 0.5),
            fractional[2] - std::floor(fractional[2] + 0.5) );
        math::rVector3d const cart(cell * centered_frac);
        math::rVector3d mini = cart;
        types::t_real dist = cart.squaredNorm();
        for(int i(-0); i < 2; ++i)
          for(int j(-0); j < 2; ++j)
            for(int k(-0); k < 2; ++k)
            {
              math::rVector3d const nvec(cart + cell * math::rVector3d(i, j, k));
              types::t_real const nnorm = nvec.squaredNorm();
              if( dist > nnorm ) 
              {
                dist = nnorm;
                mini = nvec;
              }
            }
        char const * const ptr_out = (char const * const)(PyArray_ITER_DATA(iterator.ptr()));
        *(types::t_real*)(ptr_out)           = mini[0];
        *(types::t_real*)(ptr_out + stride)  = mini[1];
        *(types::t_real*)(ptr_out + stride2) = mini[2];
      }
    }; 

    struct GaussianProj : public ToVoronoi
    {
      types::t_real sigma;
      GaussianProj  ( math::rMatrix3d const &_cell,
                      math::rVector3d const &_center, 
                      types::t_real _sigma )
                   : ToVoronoi(_cell, _center), sigma(_sigma) {}
      void init(bp::object _array)
      {
        PyArrayObject const * const ptr_array = mn::get_pyarray_pointer(_array);
        npy_intp nd = ptr_array->nd - 1;
        std::vector<npy_intp> dims(nd);
        for(size_t i(0); i < nd; ++i) dims[i] = ptr_array->dimensions[i];
        bool const isf( PyArray_ISFORTRAN(ptr_array->ob_type));
        PyObject *ptr_result = PyArray_ZEROS(nd, &(dims[0]), mn::type<types::t_real>::value, isf);
        if( ptr_result == NULL or PyErr_Occurred() != NULL ) 
        {
          if(PyErr_Occurred() == NULL)
            PyErr_SetString(PyExc_RuntimeError, "Could not create numpy array");
          boost::python::throw_error_already_set();
          return;
        }
        result = bp::object(bp::handle<>(ptr_result));
        iterator = bp::object(bp::handle<>(PyArray_IterNew(result.ptr())));
      }
      void operator()(math::rVector3d const &_orig) const
      {
        math::rVector3d const fractional(invcell * (_orig-center));
        math::rVector3d const centered_frac(
            fractional[0] - std::floor(fractional[0] + 0.5),
            fractional[1] - std::floor(fractional[1] + 0.5),
            fractional[2] - std::floor(fractional[2] + 0.5) );
        math::rVector3d const cart(cell * centered_frac);
        math::rVector3d mini = cart;
        types::t_real dist = cart.squaredNorm();
        for(int i(-0); i < 2; ++i)
          for(int j(-0); j < 2; ++j)
            for(int k(-0); k < 2; ++k)
            {
              math::rVector3d const nvec(cart + cell * math::rVector3d(i, j, k));
              types::t_real const nnorm = nvec.squaredNorm();
              if( dist > nnorm ) 
              {
                dist = nnorm;
                mini = nvec;
              }
            }
        *(types::t_real * const)(PyArray_ITER_DATA(iterator.ptr())) = std::exp(sigma * dist);
      }
    };

    template<class T_OP, class T1>
      bp::object with_op1( bp::object const &_array, T1 const &_t1)
      { T_OP op(_t1); return unloop_array_call(_array, op); }
    template<class T_OP, class T1, class T2>
      bp::object with_op2( bp::object const &_array, T1 const &_t1, T2 const &_t2 )
      { T_OP op(_t1, _t2); return unloop_array_call(_array, op); }
    template<class T_OP, class T1, class T2, class T3>
      bp::object with_op3( bp::object const &_array, T1 const &_t1, T2 const &_t2, T3 const &_t3 )
      { T_OP op(_t1, _t2, _t3); return unloop_array_call(_array, op); }




    void expose_misc()
    {
      import_array(); // needed for NumPy 
      bp::def( "to_cell", &with_op2<ToCell, math::rMatrix3d, math::rVector3d>,
               (bp::arg("positions"), bp::arg("cell"), bp::arg("center") = math::rVector3d(0,0,0) ),
               "Centers the positions in the unit cell with origin center.");
      bp::def( "to_origin", &with_op2<ToOrigin, math::rMatrix3d, math::rVector3d>,
               (bp::arg("positions"), bp::arg("cell"), bp::arg("center") = math::rVector3d(0,0,0) ),
               "Centers positions around the center in fractional coordinates.\n\n"
               "This type of centering may yield points outside the Voronoi cell of center."\
               "Use `to_voronoi` if you want the point in the Voronoi "
               "(a.k.a. Wigner-Seitz a.k.a. Brillouin) cell." );
      bp::def( "to_voronoi", &with_op2<ToVoronoi, math::rMatrix3d, math::rVector3d>,
               (bp::arg("positions"), bp::arg("cell"), bp::arg("center") = math::rVector3d(0,0,0) ),
               "Centers the positions into the voronoi "\
               "cell of center, as defined by cell. " );
      bp::def( "_gaussian_projector_impl",
               &with_op3<GaussianProj, math::rMatrix3d, math::rVector3d, types::t_real>,
               ( bp::arg("positions"), bp::arg("cell"), 
                 bp::arg("center") = math::rVector3d(0,0,0), bp::arg("alpha") = 1e0 ),
               "Computes gaussian projector around given center in unit cell.");
    }

  }
} // namespace LaDa
