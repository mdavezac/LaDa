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
    //! Calls functor _op on vector, as given by last axis.
    template<class T_OP> 
      bp::object unloop_array_call(bp::object const &_array, T_OP const &_op)
      {
        namespace mn = math::numpy;
        mn::check_is_array(_array);
        PyArrayObject *ptr_array = mn::get_pyarray_pointer(_array);
        if(PyArray_DIM(_array.ptr(), ptr_array->nd-1) != 3)
        {
          PyErr_SetString(PyExc_ValueError, "Last dimension of array should be 3.");
          bp::throw_error_already_set();
          return bp::object();
        }
        bp::object result = mn::from_template(_array, mn::type<types::t_real>::value);
  
        int last = std::max(ptr_array->nd - 1, 0);
        PyObject *i_in  = PyArray_IterAllButAxis(_array.ptr(), &last);
        PyObject *i_out = PyArray_IterAllButAxis(result.ptr(), &last);
        int const in_stride = mn::get_strides_pointer(_array)[last];
        int const in_stride2 = in_stride << 1;
        int const out_stride = mn::get_strides_pointer(result)[last];
        int const out_stride2 = out_stride << 1;
        int const type = PyArray_TYPE(_array.ptr());
        while( PyArray_ITER_NOTDONE(i_in) and PyArray_ITER_NOTDONE(i_out) )
        {
          math::rVector3d const orig( 
             mn::to_type<types::t_real>(i_in, type),
             mn::to_type<types::t_real>(i_in + in_stride, type),
             mn::to_type<types::t_real>(i_in + in_stride2, type) );
          math::rVector3d const vec( _op(orig) );
  
          char* ptr_data = (char*)(PyArray_ITER_DATA(i_out));
          *(types::t_real*)(ptr_data)               = vec[0];
          *(types::t_real*)(ptr_data + out_stride)  = vec[1];
          *(types::t_real*)(ptr_data + out_stride2) = vec[2];
          PyArray_ITER_NEXT(i_in); 
          PyArray_ITER_NEXT(i_out);
        }
        Py_DECREF(i_in);
        Py_DECREF(i_out);
        return result;
      }

    struct ToCell
    {
      math::rVector3d center;
      math::rMatrix3d cell, invcell;
      ToCell   (math::rVector3d const &_center, math::rMatrix3d const &_cell)
             : center(_center), cell(_cell) { invcell = cell.inverse(); }
      math::rVector3d operator()(math::rVector3d const &_orig) const
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
        return cell * O1;
      }
    }; 

    struct ToVoronoi : public ToCell
    {
      ToVoronoi   (math::rVector3d const &_center, math::rMatrix3d const &_cell)
                : ToCell(_center, _cell) {}
      math::rVector3d operator()(math::rVector3d const &_orig) const
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
        return mini;
      }
    }; 

    template<class T_OP>
      bp::object with_op( bp::object const &_array,
                          math::rVector3d const &_center,
                          math::rMatrix3d const &_cell )
      { return unloop_array_call(_array, T_OP(_center, _cell)); }

    void expose_misc()
    {
      bp::def( "to_cell", &with_op<ToCell>,
               (bp::arg("positions"), bp::arg("center"), bp::arg("cell")),
               "Centers the positions in the unit cell with origin center.");
      bp::def( "to_voronoi", &with_op<ToVoronoi>,
               (bp::arg("positions"), bp::arg("center"), bp::arg("cell")),
               "Centers the positions into the voronoi "\
               "cell of center, as defined by cell. " );
    }

  }
} // namespace LaDa
