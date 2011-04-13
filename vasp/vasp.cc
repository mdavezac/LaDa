#include "LaDaConfig.h"
#include "FCMangle.h"

#ifdef LADA_MPI
# include <mpi.h>
# include <dlfcn.h>

# include <boost/mpi/communicator.hpp>
#endif

#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/str.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/tuple.hpp>

namespace bp = boost::python;
#ifdef LADA_MPI
  extern "C"
  { 
    void FC_GLOBAL(vasplib, VASPLIB)(MPI_Fint *); 
    void FC_GLOBAL_(vasp_version, VASP_VERSION)(int *, int *, int *);
  }    

  # define QUOTEME_(x) #x
  # define QUOTEME(x) QUOTEME_(x)
    char const * vasp_name = QUOTEME( FC_GLOBAL(vasplib, VASPLIB) );
    char const * version_name = QUOTEME( FC_GLOBAL_(vasp_version, VASP_VERSION) );
  # undef QUOTEME 
  # undef QUOTEME_

  void *open_library(std::string const &_library)
  {
    void *dl_ptr = dlopen(_library.c_str(), RTLD_LAZY);
    if(dl_ptr == NULL)
    {
      std::ostringstream sstr;
      std::string error = dlerror();
      sstr << "Could not load escan library " << _library << ".\n"
           << "Error: " << error << "\n";
      PyErr_SetString(PyExc_ImportError, sstr.str().c_str());
      boost::python::throw_error_already_set();
      return NULL;
    }  
    return dl_ptr;
  }

  void* get_function(void *_lib, std::string const &_name)
  {
    void *call_ptr = dlsym(_lib, _name.c_str());
    if(call_ptr == NULL)
    {
      std::ostringstream sstr;
      std::string error = dlerror();
      sstr << "Could not load function " << _name << " from escan library.\n"
           << "Error: " << error << "\n";
      PyErr_SetString(PyExc_ImportError, sstr.str().c_str());
      boost::python::throw_error_already_set();
      return NULL;
    }
    return call_ptr;
  }
                                                                          
  void vasp(MPI_Comm const * const _c, std::string const &_library)                               
  {                                                                       
    void *dl_ptr = open_library(_library);
    void (*call_ptr)(MPI_Fint *) = (void(*)(MPI_Fint *)) get_function(dl_ptr, vasp_name);

    Py_BEGIN_ALLOW_THREADS;                                               
    MPI_Fint __commF = MPI_Comm_c2f( *_c );                               
    (*call_ptr)(&__commF);                        
    Py_END_ALLOW_THREADS;                                                 
    
    dlclose(dl_ptr);
  }                                                                       
                                                                          
  void boost_vasp(boost::mpi::communicator const &_c, std::string const &_library)                
  {                                                                       
    void *dl_ptr = open_library(_library);
    void (*call_ptr)(MPI_Fint *) = (void(*)(MPI_Fint *)) get_function(dl_ptr, vasp_name);

    Py_BEGIN_ALLOW_THREADS;                                               
    MPI_Comm __commC = (MPI_Comm) ( _c ) ;                                
    MPI_Fint __commF = MPI_Comm_c2f( __commC );                           
    (*call_ptr)(&__commF);                        
    Py_END_ALLOW_THREADS;                                                 

    dlclose(dl_ptr);
  }                                                                       

  bp::tuple version(std::string const &_library)
  {
    void *dl_ptr = open_library(_library);
    void (*call_ptr)(int *, int *, int *)
      = (void(*)(int *, int *, int *)) get_function(dl_ptr, version_name);

    int minor=0, major=0, medium=0;
    Py_BEGIN_ALLOW_THREADS;                                               
    (*call_ptr)(&major, &medium, &minor);                        
    Py_END_ALLOW_THREADS;                                                 

    dlclose(dl_ptr);
    return bp::make_tuple(major, medium, minor);
  }
#else
  void vasp(bp::object const &, bp::object const &)
  {
    PyErr_SetString(PyExc_RuntimeError, "LaDa compiled without MPI. Cannot run Vasp.");
    bp::throw_error_already_set();
  }
  bp::tuple version(bp::object const&) { return bp::make_tuple(5, -1, -1); } 
#endif
                                                                          
void expose_vasp()
{                                                                       
# ifdef LADA_MPI
    bp::def("vasp", &boost_vasp, (bp::arg("comm"), bp::arg("library")));           
# endif
  bp::def                                                               
  (                                                                     
    "vasp",                                                          
    &vasp,                                                         
    (bp::arg("comm"), bp::arg("library")),                                            
    "Makes a call to vasp, given the handle to a communicator.\n\n"      
    ":Parameters:\n"
    "  comm\n    The processes which vasp will use.\n"         
    "  library : str\n    A string containing the name of the vasp library.\n\n"
    "POSCAR, INCAR, IBZKPT, and other input files are expected in the current "
    "working directory.\n"
  );                                                                    
  bp::def                                                               
  (                                                                     
    "vasp",                                                          
    &vasp,                                                         
    (bp::arg("comm"), bp::arg("library")),                                            
#   ifdef LADA_MPI
      "Makes a call to vasp, given the handle to a communicator.\n\n"      
      ":Parameters:\n"
      "  comm\n    The processes which vasp will use.\n"         
      "  library : str\n    A string containing the name of the vasp library.\n\n"
      "POSCAR, INCAR, IBZKPT, and other input files are expected in the current "
      "working directory.\n"
#   else
      "LaDa compiled without mpi. Throws RuntimeError."
#   endif
  );                                                                    
  bp::def                                                               
  (                                                                     
    "version",                                                          
    &version,                                                         
    bp::arg("library"),                                            
#   ifdef LADA_MPI
      "Figures out version of vasp from library.\n\n"
      ":Parameters:\n"
      "  library : str\n    A string containing the name of the vasp library.\n"
#   else
      "LaDa compiled without MPI. Returns Fake 5.-1.-1 version."
#   endif
  );                                                                    
}

BOOST_PYTHON_MODULE(_vasp)
{
  bp::scope scope;
  scope.attr("__doc__") = "VASP-as-library extension.";
  scope.attr("__docformat__") = "restructuredtext en";
  bp::docstring_options doc_options(true, false);
  expose_vasp();
}
