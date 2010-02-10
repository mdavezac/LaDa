#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <mpi.h>
#include <boost/mpi/communicator.hpp>

#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>

extern "C" { void FC_FUNC(vasplib, VASPLIB)(MPI_Fint *); }

void vasp(MPI_Comm const * const _c) 
{
  Py_BEGIN_ALLOW_THREADS;
  MPI_Fint __commF = MPI_Comm_c2f( *_c );
  FC_FUNC(vasplib, VASPLIB)(&__commF);
  Py_END_ALLOW_THREADS;
}

void vasp2(boost::mpi::communicator const &_c) 
{ 
  Py_BEGIN_ALLOW_THREADS;
  MPI_Comm __commC = (MPI_Comm) ( _c ) ;
  MPI_Fint __commF = MPI_Comm_c2f( __commC );
  FC_FUNC(vasplib, VASPLIB)(&__commF);
  Py_END_ALLOW_THREADS;
}

BOOST_PYTHON_MODULE(_vasp)
{
  namespace bp = boost::python;
  bp::scope scope;
  scope.attr("__doc__") = "VASP-as-library extension.";
  bp::docstring_options doc_options(true, false);

  bp::def("vasp", &vasp2, bp::arg("communicator"));
  bp::def
  ( 
    "vasp",
    &FC_FUNC(vasplib, VASPLIB),
    bp::arg("communicator"),
    "Make a call to vasp, given the handle to a communicator.\n\n" 
    "@param communicator: an mpi communicator. It can be either a boost.mpi, or "
    "a mpi4py communicator."
  );
}
