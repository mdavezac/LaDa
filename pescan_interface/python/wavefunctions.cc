#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <unistd.h>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <boost/python/def.hpp>
#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/borrowed.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/errors.hpp>

#include <python/numpy_types.h>
#include <opt/initial_path.h>

#include "../interface.h"
#include "wavefunctions.hpp"


extern "C"
{
  void FC_FUNC_(escan_getwfn_datadims, ESCAN_GETWFN_DATADIMS)(int*, int*, int*);
  void FC_FUNC_(escan_copy_wfndata, ESCAN_COPY_WFNDATA)(double*, double*, int*, int*, int*);
  void FC_FUNC_(escan_read_wfns, ESCAN_COPY_WFNDATA)(int*, char const*, int*, int*, MPI_Fint*);
  void FC_FUNC_(escan_get_nr, ESCAN_GET_NR)(int *);
  void FC_FUNC_(escan_get_mr, ESCAN_GET_MR)(int *);
  void FC_FUNC_(escan_get_cell, ESCAN_GET_CELL)(double *, double *, double *);
}

namespace LaDa
{
  namespace python
  {
    void escan_get_cell( math::rMatrix3d &_cell )
    {
      double a[9], *b;
      FC_FUNC_(escan_get_cell, ESCAN_GET_CELL)(a, a+3, a+6);
      for(size_t i(0); i < 9; ++i, ++b) _cell(i%3, i/3) = *b;
    }

    namespace bp = boost::python;
    namespace bfs = boost::filesystem;
    namespace bm = boost::mpi;
 
    bp::tuple get_wavefunctions(Pescan::Interface const &_interface, bp::object const &_object)
    {
      std::vector<int> indices;
      // extract indices.
      if( bp::len(_object) == 0 )
      {
        try { indices.push_back( bp::extract<int>(_object) ); }
        catch (...) 
        {
          PyErr_SetString(PyExc_ValueError, "Second argument should integer or sequence.\n");
          bp::throw_error_already_set();
          return bp::tuple();
        }
      }
      else 
      {
        for(size_t i(0); i < bp::len(_object); ++i)
          try { indices.push_back( bp::extract<int>(_object[i]) ); }
          catch (...) 
          {
            PyErr_SetString( PyExc_ValueError, 
                             "Second argument should integer or sequence of integers.\n" );
            bp::throw_error_already_set();
            return bp::tuple();
          }
      }
      //! Makes sure indices are in 1 ... 
      int N = _interface.escan.nbstates;
      foreach(int &i, indices)
      {
        if( i < 0 ) i += _interface.escan.nbstates;
        if( i < 0 or i >= _interface.escan.nbstates)
        {
          PyErr_SetString(PyExc_IndexError, "Wavefunction index out-of-range.\n");
          bp::throw_error_already_set();
          return bp::tuple();
        }
        i += 1; // fortran index.
      }
 
      // change working directory.
      Pescan::Interface::t_Path rootpath = opt::InitialPath::path()/_interface.get_dirname();
      Pescan::Interface::t_Path origpath = bfs::current_path();
      chdir(rootpath.string().c_str());
      // reads wavefunctions
      std::string const orig =   _interface.escan.filename.string()
                               + "."
                               + boost::lexical_cast<std::string>(bm::communicator().rank());
      int a(orig.size()), b(indices.size());
      MPI_Comm __commC = (MPI_Comm) ( _interface.comm() ) ;
      MPI_Fint __commF = MPI_Comm_c2f( __commC );
      // gets dimensions.
      FC_FUNC_(escan_prepare_params, ESCAN_PREPARE_PARAMS)(&n0, &n1, &n2);
      int n0, n1, n2;
      FC_FUNC_(escan_getwfn_datadims, ESCAN_GETWFN_DATADIMS)(&n0, &n1, &n2);
      // Creates numpy objects.
      npy_intp wfn_stride[3] = { n2, n1, n0 }; // fortran dims.
      npy_intp gpoint_stride[2] = { 3, n0 };   // fortran dims
      PyObject *wfns = PyArray_SimpleNew(3, wfn_stride, NPY_CDOUBLE);
      PyObject *gpoints 
         = PyArray_SimpleNew(2, gpoint_stride, math::numpy::type<types::t_real>::value);
      char * const cptr_wfns( reinterpret_cast<PyArrayObject*>(wfns)->data );
      char * const cptr_gps( reinterpret_cast<PyArrayObject*>(gpoints)->data );
      double * const ptr_wfns( reinterpret_cast<double*>(cptr_wfns) );
      double * const ptr_gps( reinterpret_cast<double*>(cptr_gps) );
      FC_FUNC_(escan_read_wfns, ESCAN_READ_WFNS)
              ( 
                &a, orig.c_str(), 
                &b, &(indices[0]), 
                ptr_wfns, n0, n1, n2,
                ptr_gps, n0, 3,
                &__commF
              );
      
      // returns to original working directory.
      chdir(origpath.string().c_str());
      
      // returns a 2-tuple.
      return bp::make_tuple( bp::object(bp::handle<>(bp::borrowed(gpoints))), 
                             bp::object(bp::handle<>(bp::borrowed(wfns))) );
    }

    class Position
    {
      public:
        Position(bm::communicator const &_comm) : comm_(_comm) 
        {
          types::t_real data[9];
          escan_get_cell(mesh_);
          npy_intp dims[1] = {3};
          PyObject *ob = PyArray_SimpleNew(1, dims, math::numpy::type<types::t_real>::value);
          if( ob != NULL and Py_None !=  PyErr_Occurred() ) 
            reinterpret_cast<PyArrayObject*>(ob)->flags |= NPY_CARRAY_RO;
          numpy_array_ = bp::object( bp::handle<>(bp::borrowed(ob)) );
          FC_FUNC_(escan_get_nr, ESCAN_GET_NR)(&nr_);
          FC_FUNC_(escan_getwfn_datadims, ESCAN_GETWFN_DATADIMS)(&n1_, &n2_, &n3_);
          max_ = nr_ / comm_.size();
          index_ = -1;
        }
        int __len__() const { return nr_ / comm_.size(); }
        bp::object __getitem__( int _i )
        {
          if( _i < 0 ) _i += max_;
          if( _i < 0 or _i >= max_ )
          {
            PyErr_SetString(PyExc_IndexError, "index out-of-range.\n");
            bp::throw_error_already_set();
            return bp::object();
          }
          if( numpy_array_.ptr() == Py_None )
          {
            PyErr_SetString(PyExc_RuntimeError, "Could not create numpy array.\n");
            bp::throw_error_already_set();
            return bp::object();
          }
          PyArrayObject * const array = reinterpret_cast<PyArrayObject*>(numpy_array_.ptr());
          double * const data = (double*) array->data;
          int const u( (_i + 1) * comm_.rank() * nr_ / comm_.size() - 1);
          math::rVector3d const a(u/(n3_*n2_), (u%(n3_*n2_)) / n3_, (u%(n3_*n2_)) % n3_);
          math::rVector3d const b( mesh_ * a ); 
          data[0] = b(0);
          data[1] = b(1);
          data[2] = b(2);
          return numpy_array_;
        }
        void iter() const {};
        bp::object next()
        { 
          ++index_;
          if( index_ < max_ ) return __getitem__(index_);
          
          PyErr_SetString(PyExc_StopIteration, "end of range.");
          bp::throw_error_already_set();
          return bp::object();
        }

      private:
        bm::communicator comm_;
        bp::object numpy_array_;
        math::rMatrix3d mesh_;
        int index_;
        int max_;
        int nr_;
        int n3_, n2_, n1_;
    };

    bp::tuple to_realspace( bp::object const &_gwfns, bm::communicator const &_comm )
    {
      PyObject *const obj_ptr = _gwfns.ptr();
      // sanity checks
      if( not PyArray_Check(_gwfns.ptr()) ) 
      {
        PyErr_SetString(PyExc_ValueError, "Argument is not a numpy array.\n");
        bp::throw_error_already_set();
        return bp::tuple();
      }
      if(not math::numpy::is_complex(obj_ptr))
      {
        PyErr_SetString(PyExc_ValueError, "Argument is not a numpy array.\n");
        bp::throw_error_already_set();
        return bp::tuple();
      }
      PyArrayObject *array = reinterpret_cast<PyArrayObject*>(obj_ptr);
      if(array->strides[array->nd-1] != 1) 
      {
        PyErr_SetString(PyExc_ValueError, "Argument is not a contiguous numpy array. \n");
        bp::throw_error_already_set();
        return bp::tuple();
      }
      int n;
      FC_FUNC_(escan_get_mr, ESCAN_GET_MR)(&n);
      if( array->dimensions != n ) 
      {
        PyErr_SetString(PyExc_ValueError, "Unexpected array size. \n");
        bp::throw_error_already_set();
        return bp::tuple();
      }
     
      // creates output array
      std::vector<npy_intp> dims;
      for(int i(0); i < array->nd-2; ++i) dims.push_back(i);
      dims.push_back(0);
      FC_FUNC_(escan_get_nr, ESCAN_GET_NR)(&dims.back());
      if(dims.back() % 2 != 0) 
      {
        PyExc_SetString(PyExc_RuntimeError, "Wrong array dimension in pescan.");
        bp::throw_error_already_set();
        return bp::tuple();
      }
      dims.back() >>= 1;
      PyObject *result = PyArray_SimpleNew(dims.size(), &dims[0], NPY_CDOUBLE);
      if( PyErr_Occured() != Py_None )
      {
        bp::throw_error_already_set();
        return bp::tuple();
      }
      bp::object object( bp::handle<>(bp::borrowed(result)) );

      // now loops through arrays and perform fft
      if( array->nd == 1 ) 
      {
        int sign = -1;
        FC_FUNC_(d3fft_comp, D3FFT_COMP)((double*)array->data, (double*)result->data, sign);
      } 
      else
      {
        int dim = array->nd-1;
        PyArrayIterObject *real_iter 
          = reinterpret_cast<PyArrayIterObject*>(PyArray_IterAllButAxis(result, dim));
        if( not real_iter )
        {
          PyErr_SetString(PyExc_RuntimeError, "Could not iterate.\n");
          bp::throw_error_already_set();
          return bp::tuple();
        }
        PyArrayIterObject *recip_iter 
          = reinterpret_cast<PyArrayIterObject*>(PyArray_IterAllButAxis(array, dim));
        bp::object dummyA( bp::handle<>(bp::borrowed(real_iter)) );
        if( not recip_iter )
        {
          PyErr_SetString(PyExc_RuntimeError, "Could not iterate.\n");
          bp::throw_error_already_set();
          return bp::tuple();
        }
        bp::object dummyB( bp::handle<>(bp::borrowed(recip_iter)) );
        while (real_iter->index < real_iter->size and recip_iter < recip_iter->size)  
        {
          int sign = -1;
          FC_FUNC_(d3fft_comp, D3FFT_COMP)( (double*)recip_iter->dataptr,
                                            (double*)real_iter->dataptr, sign);
          PyArray_ITER_NEXT(real_iter);
          PyArray_ITER_NEXT(recip_iter);
        }
      }
      
      // finally creates real space position data.
      return bp::tuple(object, Position(_comm));
    }
 
    void expose_wfns()
    {
      import_array();
      bp::def( "wavefunctions", &get_wavefunctions, 
               (bp::arg("escan"), bp::arg("indices")), 
               "Returns wavefunctions/G-points for a given indices.\n\n"
               "@param escan: Escan functional with which calculation were performed.\n"
               "@param indices: index or indices for which to recover the wavefunctions. "
                 "Indices of wavefunctions correspond to the eigenvalues "
                 "returned by the functional during calculation.\n"
               "@result: (G-points, wavefunctions). The first item corresponds to "
                 "a 3 by x matrix with each row a G-vector. The second item is an "
                 "spin by N by x matrix holding the N wavefuntions/spinor.\n" );
    }


  } // namespace python
} // namespace LaDa
