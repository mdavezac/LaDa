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
}

namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    namespace bfs = boost::filesystem;
 
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
      Pescan::Interface::t_Path rootpath = opt::InitialPath::path();
      Pescan::Interface::t_Path origpath = bfs::current_path();
      chdir(rootpath.string().c_str());
      // reads wavefunctions
      int a(_interface.escan.filename.string().size()), b(indices.size());
      MPI_Comm __commC = (MPI_Comm) ( _interface.comm() ) ;
      MPI_Fint __commF = MPI_Comm_c2f( __commC );
      char const* fname = _interface.escan.filename.string().c_str();
      FC_FUNC_(escan_read_wfns, ESCAN_READ_WFNS)(&a, fname, &b, &(indices[0]), &__commF);
      // gets dimensions.
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
      FC_FUNC_(escan_copy_wfndata, ESCAN_COPY_WFNDATA)(ptr_wfns, ptr_gps, &n0, &n1, &n2);
      
      // returns to original working directory.
      chdir(origpath.string().c_str());
      
      // returns a 2-tuple.
      return bp::make_tuple( bp::object(bp::handle<>(bp::borrowed(gpoints))), 
                             bp::object(bp::handle<>(bp::borrowed(wfns))) );
    }
 
    void expose_wfns()
    {
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
