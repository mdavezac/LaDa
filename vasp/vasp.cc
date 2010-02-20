#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <mpi.h>
#include <boost/mpi/communicator.hpp>

#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/errors.hpp>

namespace bp = boost::python;
#ifdef LADA_EXPORT_VASP
# error LADA_EXPORT_VASP already defined
#endif
#ifdef LADA_EXPORT_NOVASP
# error LADA_EXPORT_NOVASP already defined
#endif
#define LADA_EXPORT_VASP(a)\
  extern "C" { void FC_FUNC(vasplib ## a, VASPLIB ## a)(MPI_Fint *); }   \
                                                                         \
  void vasp ## a(MPI_Comm const * const _c)                              \
  {                                                                      \
    Py_BEGIN_ALLOW_THREADS;                                              \
    MPI_Fint __commF = MPI_Comm_c2f( *_c );                              \
    FC_FUNC(vasplib ## a, VASPLIB ## a)(&__commF);                       \
    Py_END_ALLOW_THREADS;                                                \
  }                                                                      \
                                                                         \
  void boost_vasp ## a(boost::mpi::communicator const &_c)               \
  {                                                                      \
    Py_BEGIN_ALLOW_THREADS;                                              \
    MPI_Comm __commC = (MPI_Comm) ( _c ) ;                               \
    MPI_Fint __commF = MPI_Comm_c2f( __commC );                          \
    FC_FUNC(vasplib ## a, VASPLIB ## a)(&__commF);                       \
    Py_END_ALLOW_THREADS;                                                \
  }                                                                      \
                                                                         \
  void expose_vasp ## a(bp::scope &_scope)                               \
  {                                                                      \
    _scope.attr( "is_vasp" #a ) = a;                                     \
    bp::def("vasp", &boost_vasp ## a, bp::arg("communicator"));          \
    bp::def                                                              \
    (                                                                    \
      "vasp" #a,                                                         \
      &vasp ## a,                                                        \
      bp::arg("communicator"),                                           \
      "Make a call to vasp, given the handle to a communicator.\n\n"     \
      "@param communicator: the processes which vasp will use.\n"        \
      "@type communicator: boost.mpi.Communicator or mpi4py equivalent"  \
    );                                                                   \
  }

#define LADA_EXPORT_NOVASP(a) void expose_vasp ## a(bp::scope &)  {}                                 

#ifdef LADA_VASP_FOUR
  LADA_EXPORT_VASP(4)
#else
  LADA_EXPORT_NOVASP(4)
#endif
#ifdef LADA_VASP_FIVE
  LADA_EXPORT_VASP(5)
#else
  LADA_EXPORT_NOVASP(5)
#endif
#undef LADA_EXPORT_VASP

BOOST_PYTHON_MODULE(_vasp)
{
  bp::scope scope;
  scope.attr("__doc__") = "VASP-as-library extension.";
  bp::docstring_options doc_options(true, false);
  expose_vasp4(scope);
  expose_vasp5(scope);
}
