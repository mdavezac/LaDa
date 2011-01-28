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
#include "../lattice.h"
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
        math::rVector3d const vec( Crystal::into_cell(_orig-center, cell, invcell) );
        char const * const ptr_out = (char const * const)(PyArray_ITER_DATA(iterator.ptr()));
        typedef mn::type<types::t_real>::np_type t_type;
        *(t_type*)(ptr_out)           = static_cast<t_type>(vec[0]);
        *(t_type*)(ptr_out + stride)  = static_cast<t_type>(vec[1]);
        *(t_type*)(ptr_out + stride2) = static_cast<t_type>(vec[2]);
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
        math::rVector3d const vec( Crystal::zero_centered(_orig-center, cell, invcell) );
        char const * const ptr_out = (char const * const)(PyArray_ITER_DATA(iterator.ptr()));
        typedef mn::type<types::t_real>::np_type t_type;
        *(t_type*)(ptr_out)           = static_cast<t_type>(vec[0]);
        *(t_type*)(ptr_out + stride)  = static_cast<t_type>(vec[1]);
        *(t_type*)(ptr_out + stride2) = static_cast<t_type>(vec[2]);
      }
    }; 

    struct ToVoronoi : public ToCell
    {
      ToVoronoi   (math::rMatrix3d const &_cell, math::rVector3d const &_center)
                : ToCell(_cell, _center) {}
      void operator()(math::rVector3d const &_orig) const
      {
        math::rVector3d const mini( Crystal::into_voronoi(_orig-center, cell, invcell) );
        char const * const ptr_out = (char const * const)(PyArray_ITER_DATA(iterator.ptr()));
        typedef mn::type<types::t_real>::np_type t_type;
        *(t_type*)(ptr_out)           = static_cast<t_type>(mini[0]);
        *(t_type*)(ptr_out + stride)  = static_cast<t_type>(mini[1]);
        *(t_type*)(ptr_out + stride2) = static_cast<t_type>(mini[2]);
      }
    }; 

    //! Creates gaussian projector in r-space.
    struct GaussianProj : public ToCell
    {
      types::t_real sigma;
      GaussianProj  ( math::rMatrix3d const &_cell,
                      math::rVector3d const &_center, 
                      types::t_real _sigma )
                   : ToCell(_cell, _center), sigma(_sigma) {}
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
        math::rVector3d const mini( Crystal::into_voronoi(_orig-center, cell, invcell) );
        types::t_real const dist = mini.squaredNorm();
        typedef mn::type<types::t_real>::np_type t_type;
        *(t_type*)(PyArray_ITER_DATA(iterator.ptr())) = static_cast<t_type>(std::exp(sigma * dist));
      }
    };


    //! Returns boolean arrays indicating points on lattice.
    struct IsOnLattice : public ToCell
    {
      types::t_real tolerance;
      IsOnLattice  ( math::rMatrix3d const &_cell,
                     math::rVector3d const &_center, 
                     types::t_real _tolerance )
                  : ToCell(_cell, _center), tolerance(_tolerance) {}
      void init(bp::object _array)
      {
        PyArrayObject const * const ptr_array = mn::get_pyarray_pointer(_array);
        npy_intp nd = ptr_array->nd - 1;
        std::vector<npy_intp> dims(nd);
        for(size_t i(0); i < nd; ++i) dims[i] = ptr_array->dimensions[i];
        bool const isf( PyArray_ISFORTRAN(ptr_array->ob_type));
        PyObject *ptr_result = PyArray_ZEROS(nd, &(dims[0]), mn::type<bool>::value, isf);
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
        math::rVector3d const vec(Crystal::into_cell(_orig-center, cell, invcell));
        typedef mn::type<bool>::np_type t_type;
        *(t_type*)(PyArray_ITER_DATA(iterator.ptr()))
            = static_cast<t_type>(vec.squaredNorm() < tolerance);
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
      bp::def( "is_on_lattice",
               &with_op3<IsOnLattice, math::rMatrix3d, math::rVector3d, types::t_real>,
               ( bp::arg("positions"), bp::arg("cell"), 
                 bp::arg("origin") = math::rVector3d(0,0,0),
                 bp::arg("tolerance") = types::tolerance ),
               "Computes whether a set of points is on a lattice defined by origin and cell." );
    }

  }
} // namespace LaDa
