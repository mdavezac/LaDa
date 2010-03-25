#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <unistd.h>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <boost/python/def.hpp>
#include <boost/python/object.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/borrowed.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/return_arg.hpp>

#include <python/numpy_types.h>
#include <opt/initial_path.h>

#include "../interface.h"
#include "wavefunctions.hpp"


extern "C"
{
  void FC_FUNC_(escan_wfns_init, ESCAN_WFNS_INIT)(int*, char const*, MPI_Fint*);
  void FC_FUNC_(escan_wfns_get_array_dimensions, ESCAN_WFNS_GET_ARRAY_DIMENSIONS)(int*, int*, int*);
  void FC_FUNC_(escan_wfns_read, ESCAN_WFNS_READ)(int*, int*, int*, int*, int const*, double*, double*);
  void FC_FUNC_(escan_get_nr, ESCAN_GET_NR)(int *);
  void FC_FUNC_(escan_get_mr_n, ESCAN_GET_MR_N)(int *);
  void FC_FUNC_(escan_get_n2_n3, ESCAN_GET_N2_N3)(int*, int*);
  void FC_FUNC_(escan_get_cell, ESCAN_GET_CELL)(double *, double *, double *);
  void FC_FUNC_(d3fft_comp, D3FFT_COMP)(double*, double*, int*);
  void FC_FUNC_(escan_wfns_cleanup, ESCAN_WFNS_CLEANUP)();
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
 
    class Wfns
    {
      public:
        Wfns(Pescan::Interface const &_pescan, bp::object const &_object);
         
        bp::tuple __enter__() const;
        void __exit__(bp::object const&, bp::object const&, bp::object const&) const 
          { FC_FUNC_(escan_wfns_cleanup, ESCAN_WFNS_CLEANUP)(); }

      private:
        Pescan::Interface interface_;
        std::vector<int> indices_;
    };

    Wfns :: Wfns   (Pescan::Interface const &_pescan, bp::object const &_object) 
                 : interface_(_pescan)
    {
      // extract indices.
      if( bp::len(_object) == 0 )
      {
        try { indices_.push_back( bp::extract<int>(_object) ); }
        catch (...) 
        {
          PyErr_SetString(PyExc_ValueError, "Second argument should integer or sequence.\n");
          bp::throw_error_already_set();
          return;
        }
      }
      else 
      {
        for(size_t i(0); i < bp::len(_object); ++i)
          try { indices_.push_back( bp::extract<int>(_object[i]) ); }
          catch (...) 
          {
            PyErr_SetString( PyExc_ValueError, 
                             "Second argument should integer or sequence of integers.\n" );
            bp::throw_error_already_set();
            return;
          }
      }
      //! Makes sure indices are in 1 ... 
      int N = interface_.escan.nbstates;
      foreach(int &i, indices_)
      {
        if( i < 0 ) i += interface_.escan.nbstates;
        if( i < 0 or i >= interface_.escan.nbstates)
        {
          PyErr_SetString(PyExc_IndexError, "Wavefunction index out-of-range.\n");
          bp::throw_error_already_set();
          return;
        }
        i += 1; // fortran index.
      }
    } 
    bp::tuple Wfns::__enter__() const 
    {
      // change working directory.
      Pescan::Interface::t_Path rootpath = opt::InitialPath::path()/interface_.get_dirname();
      Pescan::Interface::t_Path origpath = bfs::current_path();
      chdir(rootpath.string().c_str());
      // prepares to read wavefinctions
      std::string const orig =   interface_.escan.filename.string()
                               + "."
                               + boost::lexical_cast<std::string>(bm::communicator().rank());
      int a(orig.size()), b(indices_.size());
      MPI_Comm __commC = (MPI_Comm) ( interface_.comm() ) ;
      MPI_Fint __commF = MPI_Comm_c2f( __commC );
      FC_FUNC_(escan_wfns_init, ESCAN_WFNS_INIT)(&a, orig.c_str(), &__commF);
      // gets dimensions.
      int n0, n1(indices_.size()), n2, g0, g1(3);
      FC_FUNC_(escan_wfns_get_array_dimensions, ESCAN_WFNS_GET_ARRAY_dimensions)(&n0, &n2, &g0);
      
      // Creates numpy objects.
      npy_intp wfn_dims[3] = { n0, n1, n2 }; // fortran dims.
      npy_intp gpoint_dims[2] = { g0, g1 };   // fortran dims
      PyObject *wfns = PyArray_ZEROS(3, wfn_dims, NPY_CDOUBLE, 1);
      PyObject *gpoints 
         = PyArray_ZEROS(2, gpoint_dims, math::numpy::type<types::t_real>::value, 1);
      char * const cptr_wfns( reinterpret_cast<PyArrayObject*>(wfns)->data );
      char * const cptr_gps( reinterpret_cast<PyArrayObject*>(gpoints)->data );
      double * const ptr_wfns( reinterpret_cast<double*>(cptr_wfns) );
      double * const ptr_gps( reinterpret_cast<double*>(cptr_gps) );
      
      reinterpret_cast<PyArrayObject*>(wfns)->dimensions[0]
        = reinterpret_cast<PyArrayObject*>(gpoints)->dimensions[0];

      // finally reads wavefunctions
      FC_FUNC_(escan_wfns_read, ESCAN_WFNS_READ)
              ( 
                &n0, &n1, &n2, &g0,  // dimensions
                &(indices_[0]),       // indices_ to wavefunctions
                ptr_wfns,            // pointer to wavefunctions data
                ptr_gps              // pointer to gpoint  data
              );
      

      // returns to original working directory.
      chdir(origpath.string().c_str());
      
      
      // returns a 2-tuple.
      return bp::make_tuple( bp::object(bp::handle<>(wfns)), 
                             bp::object(bp::handle<>(gpoints)) );
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
          if( ob != NULL and PyErr_Occurred() != NULL ) 
            reinterpret_cast<PyArrayObject*>(ob)->flags |= NPY_CARRAY_RO;
          numpy_array_ = bp::object( bp::handle<>(ob) );
          FC_FUNC_(escan_get_nr, ESCAN_GET_NR)(&nr_);
          FC_FUNC_(escan_get_n2_n3, ESCAN_GET_N2_n3)(&n2_, &n3_);
          max_ = nr_ / comm_.size();
          index_ = -1;
        }
        int __len__() const { return max_; }
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
        int n3_, n2_;
    };

    bp::tuple to_realspace( bp::object const &_gwfns, bm::communicator const &_comm )
    {
      PyObject *const obj_ptr = _gwfns.ptr();
      // sanity checks
      if( not PyArray_Check(obj_ptr) ) 
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
      if(array->strides[0] != 16) 
      {
        PyErr_SetString(PyExc_ValueError, "Argument is not a contiguous numpy array of complexes. \n");
        bp::throw_error_already_set();
        return bp::tuple();
      }
      int n0, n2, g0;
      FC_FUNC_(escan_wfns_get_array_dimensions, ESCAN_WFNS_GET_ARRAY_DIMENSIONS)(&n0, &n2, &g0);
      if( array->nd == 1 )
      {
        if(array->dimensions[0] != n0)
        {
          PyErr_SetString(PyExc_ValueError, "Unexpected array size. \n");
          bp::throw_error_already_set();
          return bp::tuple();
        }
      }
      else if( array->strides[0] * n0 < array->strides[1] ) 
      {
        PyErr_SetString(PyExc_ValueError, "Unexpected array size along first dimension. \n");
        bp::throw_error_already_set();
        return bp::tuple();
      }
     
      // creates output array
      FC_FUNC_(escan_get_mr_n, ESCAN_GET_MR_N)(&n0);
      if(n0 % 2 != 0) 
      {
        PyErr_SetString(PyExc_RuntimeError, "Wrong array dimension in pescan.");
        bp::throw_error_already_set();
        return bp::tuple();
      }
      std::vector<npy_intp> dims(1, n0>>1);
      for(int i(1); i < array->nd; ++i) dims.push_back(array->dimensions[i]);
      PyObject *result = PyArray_ZEROS(dims.size(), &dims[0], NPY_CDOUBLE, 1);
      if( PyErr_Occurred() != NULL )
      {
        bp::throw_error_already_set();
        return bp::tuple();
      }
      bp::object resob = bp::object(bp::handle<>(result));

      // now loops through arrays and perform fft
      if( array->nd == 1 ) 
      {
        int sign = -1;
        FC_FUNC_(d3fft_comp, D3FFT_COMP)
        (
          (double*)array->data, 
          (double*)reinterpret_cast<PyArrayObject*>(result)->data,
          &sign
        );
      } 
      else
      {
        int dim = 0;
        PyArrayIterObject *real_iter 
          = reinterpret_cast<PyArrayIterObject*>(PyArray_IterAllButAxis(result, &dim));
        if( not real_iter )
        {
          PyErr_SetString(PyExc_RuntimeError, "Could not iterate.\n");
          bp::throw_error_already_set();
          return bp::tuple();
        }
        PyArrayIterObject *recip_iter 
          = reinterpret_cast<PyArrayIterObject*>
            (
              PyArray_IterAllButAxis(reinterpret_cast<PyObject*>(array), &dim)
            );
        // object should release on destruction
        bp::object dummyA = bp::object(bp::handle<>(bp::borrowed(reinterpret_cast<PyObject*>(real_iter))));
        if( not recip_iter )
        {
          PyErr_SetString(PyExc_RuntimeError, "Could not iterate.\n");
          bp::throw_error_already_set();
          return bp::tuple();
        }
        // object should release on destruction
        bp::object dummyB = bp::object(bp::handle<>(bp::borrowed(reinterpret_cast<PyObject*>(recip_iter))));
        while (real_iter->index < real_iter->size and recip_iter->index < recip_iter->size)  
        {
          int sign = -1;
          FC_FUNC_(d3fft_comp, D3FFT_COMP)( (double*)recip_iter->dataptr,
                                            (double*)real_iter->dataptr, &sign);
          PyArray_ITER_NEXT(real_iter);
          PyArray_ITER_NEXT(recip_iter);
        }
      }
      
      int nr;
      FC_FUNC_(escan_get_nr, ESCAN_GET_NR)(&nr);
      reinterpret_cast<PyArrayObject*>(result)->dimensions[0] = nr / _comm.size();
      // finally creates real space position data.
      return bp::make_tuple(resob, Position(_comm));
    }
 
    void expose_wfns()
    {
      import_array();
      bp::class_<Wfns>
      (
        "Wavefunctions", 
        "Context with temporary arrays to wavefunctions and corresponding g-vectors.\n\n"
        "@param escan: Escan functional with which calculation were performed.\n"
        "@param indices: index or indices for which to recover the wavefunctions. "
          "Indices of wavefunctions correspond to the eigenvalues "
          "returned by the functional during calculation.\n"
        "@return: (wavefunctions, g-vectors). The second item corresponds to "
          "a 3 by x matrix with each row a G-vector. The first item is an "
          "spin by N by x matrix holding the N wavefuntions/spinor.\n",
        bp::init<Pescan::Interface const&, bp::object const>()
      ).def("__enter__", &Wfns::__enter__)
       .def("__exit__", &Wfns::__exit__);
      bp::class_<Position>("Position", "Real-space wavefunction vectors.", bp::no_init)
        .def("__getitem__", &Position::__getitem__)
        .def("__len__", &Position::__len__)
        .def("__iter__", &Position::iter, bp::return_self<>())
        .def("next", &Position::next);
      bp::def( "to_realspace", &to_realspace, (bp::arg("wfns"), bp::arg("comm")),
               "Returns wavefunctions in real-space.\n\n"
               "@param wfns: Array of reciprocal-space wavefunctions.\n"
               "@type wfns: numpy array\n"
               "@param mpicomm: mpi communicator.\n"
               "@type mpicomm: boost\n"
               "@return: (numpy array real-space wfns, generator/array of positions\n");
    }


  } // namespace python
} // namespace LaDa
