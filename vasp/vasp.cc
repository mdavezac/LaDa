#include "LaDaConfig.h"

#include <mpi.h>
#include <boost/mpi/communicator.hpp>

#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/str.hpp>
#include <boost/python/errors.hpp>

namespace bp = boost::python;
extern "C"
{ 
  void FC_FUNC(vasplib, VASPLIB)(MPI_Fint *); 
  void FC_FUNC_(vasp_version, VASP_VERSION)(int *, int *, int *);
}    
                                                                        
void vasp(MPI_Comm const * const _c)                               
{                                                                       
  Py_BEGIN_ALLOW_THREADS;                                               
  MPI_Fint __commF = MPI_Comm_c2f( *_c );                               
  FC_FUNC(vasplib, VASPLIB)(&__commF);                        
  Py_END_ALLOW_THREADS;                                                 
}                                                                       
                                                                        
void boost_vasp(boost::mpi::communicator const &_c)                
{                                                                       
  Py_BEGIN_ALLOW_THREADS;                                               
  MPI_Comm __commC = (MPI_Comm) ( _c ) ;                                
  MPI_Fint __commF = MPI_Comm_c2f( __commC );                           
  FC_FUNC(vasplib, VASPLIB)(&__commF);                        
  Py_END_ALLOW_THREADS;                                                 
}                                                                       
                                                                        
void expose_vasp()
{                                                                       
  bp::def("vasp", &boost_vasp, bp::arg("communicator"));           
  bp::def                                                               
  (                                                                     
    "vasp",                                                          
    &vasp,                                                         
    bp::arg("communicator"),                                            
    "Make a call to vasp, given the handle to a communicator.\n\n"      
    "POSCAR, INCAR, IBZKPT, and other input files are expected in the current "
    "working directory.\n"
    "@param communicator: the processes which vasp will use.\n"         
    "@type communicator: boost.mpi.Communicator or mpi4py equivalent"   
  );                                                                    
}

BOOST_PYTHON_MODULE(_vasp)
{
  bp::scope scope;
  scope.attr("__doc__") = "VASP-as-library extension.";
  bp::docstring_options doc_options(true, false);
  int minor=0, major=0, medium=0;
  FC_FUNC_(vasp_version, VASP_VERSION)(&major, &medium, &minor);
  std::ostringstream sstr; sstr << major << "." << medium << "." << minor;
  scope.attr( "version" ) = sstr.str();
  expose_vasp();
}
