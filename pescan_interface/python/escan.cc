//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/class.hpp>
#include <boost/python/object.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/def.hpp>
#include <boost/python/str.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/data_members.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <opt/initial_path.h>

#include <python/numpy_types.h>
#include <vff/python/vff.hpp>
#include <opt/tuple_serialize.h>

#include "../interface.h"
#include "escan.hpp"

//! \cond
extern "C"
{
  void FC_FUNC_(iaga_set_mpi, IAGA_SET_MPI)( MPI_Fint * );
  void FC_FUNC_(getvlarg, GETVLARG)();
  void FC_FUNC_(iaga_just_call_escan, IAGA_JUST_CALL_ESCAN)();
}
//! \endcond

void just_call_escan(boost::mpi::communicator const &_c)
{
  MPI_Comm __commC = (MPI_Comm) ( _c ) ;
  MPI_Fint __commF = MPI_Comm_c2f( __commC );
  FC_FUNC_(iaga_set_mpi, IAGA_SET_MPI)( &__commF );

  FC_FUNC_(getvlarg, GETVLARG)();
  FC_FUNC_(iaga_just_call_escan, IAGA_just_CALL_ESCAN)();
}
void just_call_escan2(MPI_Comm &_comm)
{
  MPI_Fint __commF = MPI_Comm_c2f( _comm );
  FC_FUNC_(iaga_set_mpi, IAGA_SET_MPI)( &__commF );

  FC_FUNC_(getvlarg, GETVLARG)();
  FC_FUNC_(iaga_just_call_escan, IAGA_just_CALL_ESCAN)();
}

namespace LaDa
{
  namespace bp = boost::python;
  namespace python
  {
    template<class T> 
      struct pickle_escan : boost::python::pickle_suite
      {
        static  bp::tuple getinitargs(const T&) { return bp::tuple(); }

        static bp::tuple getstate(T const &_escan)
        {
          std::ostringstream ss;
          boost::archive::text_oarchive oa( ss );
          oa << _escan;

          return bp::make_tuple( ss.str() );
        }
        static void setstate(T &_escan, bp::tuple _state)
        {
          if( bp::len( _state ) != 1 )
          {
            PyErr_SetObject(PyExc_ValueError,
                            ("expected 1-item tuple in call to __setstate__; got %s"
                             % _state).ptr()
                );
            bp::throw_error_already_set();
          }
          const std::string str = bp::extract< std::string >( _state[0] );
          std::istringstream ss( str.c_str() );
          boost::archive::text_iarchive ia( ss );
          ia >> _escan;
        }
      };



    template<class T_TYPE>
      void Escan_from_XML(T_TYPE &_type, const std::string &_filename )
      {
        if( not opt::InitialPath::is_initialized() ) opt::InitialPath::init();

        TiXmlDocument doc( _filename ); 
        TiXmlHandle docHandle( &doc ); 
      
        if( not doc.LoadFile() )
        {
          PyErr_SetString
          (
            PyExc_RuntimeError,
            (   "Could not parse xml file " + _filename + ".\n" 
              + doc.ErrorDesc() + "\n" ).c_str()
          );
          bp::throw_error_already_set();
        }
        else if( not docHandle.FirstChild("Job").Element() )
        {
          PyErr_SetString
          (
            PyExc_RuntimeError,
            ("Could not find <Job> in xml file " + _filename + ".\n").c_str()
          );
          bp::throw_error_already_set();
        }
        else
        {
          try
          {
            if(not _type.Load( *docHandle.FirstChild("Job").Element() ) )
            {
              PyErr_SetString
              (
                PyExc_RuntimeError,
                ("Could not load VFF functional from xml file " + _filename + ".\n").c_str()
              );
              bp::throw_error_already_set();
            }
          }
          catch(std::exception &_e)
          {
            PyErr_SetString
            (
              PyExc_RuntimeError,
              (
                   "Encountered error while loading  VFF functional from xml file "
                 + _filename + ".\n"
                 + _e.what() + "\n"
              ).c_str()
            );
            bp::throw_error_already_set();
          }
        }
      }
    template<class T> T* cload(std::string const &_filename) 
    {
      T* result(new T);
      if( not result )
      {
        PyErr_SetString(PyExc_RuntimeError, "Memory error.\n");
        bp::throw_error_already_set();
        return NULL;
      }
      Escan_from_XML(*result, _filename); 
      if( PyErr_Occurred() != NULL ) 
      {
        delete result;
        return NULL;
      }
      return result;
    }
#   ifdef _MPI
      template<class T> T* mpi_create(boost::mpi::communicator *_c)
      {
        T* result(new T);
        if( not result )
        {
          PyErr_SetString(PyExc_RuntimeError, "Memory error.\n");
          bp::throw_error_already_set();
          return NULL;
        }
        result->set_mpi(_c);
        return result;
      }
      template<class T> T* mpi_create2( std::string const &_filename,
                                        boost::mpi::communicator *_c )
                                       
      {
        T* result(new T);
        if( not result )
        {
          PyErr_SetString(PyExc_RuntimeError, "Memory error.\n");
          bp::throw_error_already_set();
          return NULL;
        }
        result->set_mpi(_c);
        Escan_from_XML(*result, _filename); 
        if( PyErr_Occurred() != NULL ) 
        {
          delete result;
          return NULL;
        }
        return result;
      }
#   else 
      template<class T> T* mpi_create(bp::object const &)
      {
        T* result(new T);
        if( not result )
        {
          PyErr_SetString(PyExc_RuntimeError, "Memory error.\n");
          bp::throw_error_already_set();
          return NULL;
        }
        return result;
      }
      template<class T> T* mpi_create2(bp::object const &, std::string const &_filename)
        { return cload(_filename); }
#   endif
   
#   ifdef LADA_PARAMS
#     error LADA_PARAMS already exists.
#   endif
#   define LADA_PARAMS(type, name) \
      type get_ ## name (Pescan::Interface const &_interface) \
        { return _interface.escan.name; } \
      void set_ ## name(Pescan::Interface &_interface, type const &_value) \
        { _interface.escan.name = _value; }
      LADA_PARAMS(types::t_int, nbstates)
      LADA_PARAMS(types::t_int, nlines)
      LADA_PARAMS(types::t_real, tolerance)
      LADA_PARAMS(types::t_real, rcut)
      LADA_PARAMS(math::rVector3d, kpoint)
      LADA_PARAMS(Pescan::Interface::t_method, method)
      LADA_PARAMS(Pescan::Interface::Escan::t_potential, potential)
      LADA_PARAMS(types::t_real, kinscal)
#   undef LADA_PARAMS
#   define LADA_PARAMS(type, name, var) \
      type get_ ## name(Pescan::Interface const &_interface) \
        { return _interface.escan.var; } \
      void set_ ## name(Pescan::Interface &_interface, type const &_value) \
        { _interface.escan.var = _value; }
      LADA_PARAMS(types::t_real, reference, Eref)
      LADA_PARAMS(types::t_real, smoothness, smooth)
#   undef LADA_PARAMS
#   define LADA_PARAMS(_A_, _B_) \
      std::string get_ ## _B_ (LaDa::Pescan::Interface const &_a) {return _a._A_.string();} \
      void set_ ## _B_ (LaDa::Pescan::Interface &_a, const std::string &_b)  {_a._A_ = _b;} 
      LADA_PARAMS(escan.filename, input_filename)
      LADA_PARAMS(escan.output, output_filename)
      LADA_PARAMS(atom_input, vff_inputfile)
      LADA_PARAMS(stdout_file, stdout)
      LADA_PARAMS(stderr_file, stderr)
#   undef LADA_PARAMS
   
    typedef LaDa::Pescan::Interface::GenPot::t_MeshTuple t_MeshTuple;
    
    bp::tuple get_impl( t_MeshTuple const& _mesh )
    {
      namespace bt = boost::tuples;
      return bp::make_tuple( bt::get<0>(_mesh), bt::get<1>(_mesh), bt::get<2>(_mesh) );
    }
    void set_impl( bp::tuple const &_tuple, t_MeshTuple & _mesh )
    {
      namespace bt = boost::tuples;
      try
      {
        if( bp::len(_tuple) != 3 )
        {
          PyErr_SetString( PyExc_RuntimeError, "Should be a 3-tuple.\n" );
          bp::throw_error_already_set();
          return;
        }
        bt::get<0>(_mesh) = bp::extract<types::t_int>( _tuple[0] );
        bt::get<1>(_mesh) = bp::extract<types::t_int>( _tuple[1] );
        bt::get<2>(_mesh) = bp::extract<types::t_int>( _tuple[2] );
      }
      catch(...)
      {
        PyErr_SetString( PyExc_RuntimeError, "Could not translate tuple.\n" );
        bp::throw_error_already_set();
      }
    }
    
    bp::tuple get_mesh( LaDa::Pescan::Interface::GenPot const &_genpot )
      { return get_impl(_genpot.mesh); }
    void set_mesh( LaDa::Pescan::Interface::GenPot &_genpot, bp::tuple const &_t )
    {
      namespace bt = boost::tuples;
      bool const change_multcell
      (
            bt::get<0>( _genpot.mesh) == bt::get<0>(_genpot.multiple_cell)
        and bt::get<1>( _genpot.mesh) == bt::get<1>(_genpot.multiple_cell)
        and bt::get<2>( _genpot.mesh) == bt::get<2>(_genpot.multiple_cell)
      );
      set_impl(_t, _genpot.mesh); 
      if( change_multcell ) _genpot.multiple_cell = _genpot.mesh;
    }
    bp::tuple get_mcell( LaDa::Pescan::Interface::GenPot const &_genpot )
      { return get_impl(_genpot.multiple_cell); }
    void set_mcell( LaDa::Pescan::Interface::GenPot &_genpot, bp::tuple const &_t )
      { set_impl(_t, _genpot.multiple_cell); }
    bp::tuple get_sbox( LaDa::Pescan::Interface::GenPot const &_genpot )
      { return get_impl(_genpot.small_box); }
    void set_sbox( LaDa::Pescan::Interface::GenPot &_genpot, bp::tuple const &_t )
      { set_impl(_t, _genpot.small_box); }
   
    bp::str get_dir( LaDa::Pescan::Interface const &_self )
      { return bp::str( _self.get_dirname().string() ); }
    void set_dir( LaDa::Pescan::Interface &_self, bp::str const &_str )
    {
      std::string const string = bp::extract<std::string>(_str);
      return _self.set_dirname(string); 
    }
   
    bp::object get_eigenvalues(Pescan::Interface const& _interface)
    {
      namespace numpy = LaDa::math::numpy;
      npy_intp dims[1] = {_interface.eigenvalues.size()};

      PyObject *eigs = PyArray_SimpleNewFromData( 1, dims, numpy::type<types::t_real>::value,
                                                  (void*)&(_interface.eigenvalues[0]) );
      return bp::object( bp::handle<>(bp::borrowed(eigs)) );
    }
    bp::object __call__(Pescan::Interface &_interface) 
    {
      if( not _interface() )
      {
        PyErr_SetString(PyExc_RuntimeError, "Something went wrong in escan.\n");
        bp::throw_error_already_set();
      }
      return get_eigenvalues(_interface);
    }
    void set_scale(Pescan::Interface &_interface, Crystal::TStructure<std::string> const &_str)
      {  _interface.escan.scale = _str.scale; }
    template<class T>
      bp::object __call2__(Pescan::Interface &_interface,
                           T &_vff, Crystal::TStructure<std::string> const &_str)
      {
        Pescan::Interface::t_Path const path = _interface.atom_input;
        Pescan::Interface::t_Path rootdir = path.root_path();
        Pescan::Interface::t_Path filename = path.filename();
        if( filename.empty() ) filename = "atom_input";
        boost::mpi::communicator world;
        std::ostringstream sstr; sstr << world.rank();
        _interface.atom_input = rootdir / filename.replace_extension(sstr.str());
        print_escan_input(_vff, _interface.atom_input.string(), _str);
        set_scale(_interface, _str);
        if( not _interface() )
        {
          PyErr_SetString(PyExc_RuntimeError, "Something went wrong in escan.\n");
          bp::throw_error_already_set();
        }
        _interface.atom_input = path;
        return get_eigenvalues(_interface);
      }

    template<class T> size_t nb_valence_states( T const &_str ) 
      { return Pescan::nb_valence_states( _str ); }


    void expose_escan()
    {
      typedef LaDa::Pescan::Interface t_Escan;
      bp::enum_<t_Escan::t_method>( "method", "Diagonalisation method." )
        .value( "folded", t_Escan::FOLDED_SPECTRUM )
        .value( "full_diagonalization", t_Escan::ALL_ELECTRON )
        .export_values();
      bp::enum_<t_Escan::Escan::t_potential>( "potential", "Type of the hamiltonian" )
        .value( "local", t_Escan::Escan::LOCAL )
        .value( "nonlocal", t_Escan::Escan::NONLOCAL )
        .value( "spinorbit", t_Escan::Escan::SPINORBIT )
        .export_values();

      bp::class_< t_Escan >
      ( 
        "Escan", 
        "Wrapper around the nanopse-escan functional.\n\n"
        "Initialization can take:\n"
         "  - No argument.\n"
         "  - A single string argument representing the path to an XML input file.\n" 
         "  - A single boost.mpi communicator.\n"
         "  - A string argument (see above), followed by a boost.mpi communicator.\n\n" 
         "If compiled without mpi, including the communicator will have no effect.\n"
      ) .def( bp::init< t_Escan const& >() )
        .def( "__init__", bp::make_constructor(&cload<t_Escan>) )\
        .def( "__init__", bp::make_constructor(&mpi_create<t_Escan>) )\
        .def( "__init__", bp::make_constructor(&mpi_create2<t_Escan>) )\
        .add_property
        (
           "directory", &get_dir, &set_dir, 
           "Directory where calculations are performed.\n"
        )
        .add_property
        ( 
          "scale", &t_Escan::get_scale, &set_scale,
          "Internal escan scale.\n\n"
          "Prior to calculation, set as C{escan.scale = structure} where \"structure\" is a "
          "L{Structure<lada.crystal.Structure>}.\n"
        )
#       define LADA_PARAMS(name, docstring)\
          .add_property(#name, &get ## name, &set ## name, docstring)
#       undef LADA_PARAMS
#       define LADA_PARAMS(name, docstring) \
          .add_property(#name, &get_ ## name, &set_ ## name, docstring)
          LADA_PARAMS(input_filename, "Input filename for pescan.")
          LADA_PARAMS(output_filename, "Output filename for pescan.")
          LADA_PARAMS(rcut, "Real-space cutoff.")
          LADA_PARAMS(smoothness, "Smoothness factor of plane-wave cutoff.")
          LADA_PARAMS(kinscal, "Kinetic scaling parameter.")
          LADA_PARAMS(vff_inputfile, "Structure input file.")
          LADA_PARAMS(stdout, "Filename of the standard output.\n\n"
                      "If empty, then stdout is not redirected.")
          LADA_PARAMS(stderr, "Filename of the standard error.\n\n"
                      "If empty, then stderr is not redirected.")
          LADA_PARAMS(nbstates, "Number of states to compute.")
          LADA_PARAMS(method, "Diagonalization method: L{folded} or L{full_diagonalization}")
          LADA_PARAMS(potential, "Hamiltonian method: L{local}, L{spinorbit}, or L{nonlocal}.")
          LADA_PARAMS(reference, "Reference energy for folded spectrum method." )
          LADA_PARAMS(nlines, "Maximum number of line optimizations in conjugate-gradient.")
          LADA_PARAMS(tolerance, "Tolerance of the diagonalisation procedure.")
#       undef LADA_PARAMS
        .add_property
        (
          "kpoint", 
           bp::make_function(&get_kpoint, bp::return_value_policy<bp::return_by_value>()),
           &set_kpoint, "K-point, in units of M{2S{pi}/a}, at which to perform calculation."
        )
        .def_readwrite( "destroy_directory", &t_Escan::do_destroy_dir,
                        "If true, directory where calculations are carried out is destroyed "
                        "at the end of a calculation." )
        .add_property
        ( 
          "eigenvalues", 
          bp::make_function(&get_eigenvalues, bp::with_custodian_and_ward_postcall<1,0>()),
          "Computed eigenvalues as a read-only numpy array of real values." 
        )
        .def_readwrite("genpot", &t_Escan::genpot, "(L{GenPot}) Parameters of the atomic "
            "potentials.\n")
        .def_readwrite("verbose", &t_Escan::verbose, "Verbose pescan output on true.")
        .def( "fromXML",  &Escan_from_XML<t_Escan>, bp::arg("file"),
              "Loads escan parameters from an XML file." )
        .def
        (
          "__call__", &__call2__<t_Vff>, 
          (bp::arg("vff"), bp::arg("structure")),
          bp::with_custodian_and_ward_postcall<1,0>() 
        )
        .def
        ( 
          "__call__", &__call2__<t_LayeredVff>,
          (bp::arg("vff"), bp::arg("structure")),
          "Performs a calculation.\n\n" 
          "The usual sequence of call is as follows:\n\n"
          ">>> # First initialize escan and vff"
          ">>> escan = lada.escan.Escan(\"input.xml\", boost.mpi.world)\n"
          ">>> vff = lada.vff.Vff(\"input.xml\", boost.mpi.world)\n"
          ">>> # Finally, calls escan functional. \n"
          ">>> # Other parameters (reference, kpoint...) can be changed prior to call.\n"
          ">>> eigenvalues = escan(vff, structure)\n\n"
          "Escan needs to be passed a valid vff structure to compute the microstrain.\n\n"
          "@param vff: A vff functional.\n"
          "@type vff: L{vff.Vff} or L{vff.LayeredVff}\n"
          "@param structure: The structure for which to perform calculations.\n"
          "@type structure: L{crystal.Structure}\n"
          "@return: numpy vector also available as self.L{eigenvalues}. As "
          "such, it will change from one calculation to the other. Copy if you "
          "want to store the results.\n",
          bp::with_custodian_and_ward_postcall<1,0>() 
        )
        .def( "set_mpi", &t_Escan::set_mpi, "Sets the boost.mpi communicator." )
        .def_pickle( pickle_escan< t_Escan >() );

      bp::def("_call_escan", &just_call_escan);
      bp::def("_call_escan2", &just_call_escan, "Private interface. @see lada.escan.call_escan.");

      bp::def( "nb_valence_states", &nb_valence_states<Crystal::TStructure<std::string> >,
               bp::arg("structure"), "Returns the number of valence states in a structure." );
    }

    void expose_genpot()
    {

      typedef LaDa::Pescan::Interface::GenPot t_GenPot;
      import_array();
      bp::class_< t_GenPot >( "GenPot", "Wrapper around genpot parameters.", bp::no_init ) 
        .add_property( "mesh", &get_mesh, &set_mesh, "Tuple defining the FFT mesh.\n" ) 
        .add_property( "multiple_cell", &get_mcell, &set_mcell,
                       "Tuple defining the cell-grid for large computations.\n" ) 
        .add_property( "small_box", &get_sbox, &set_sbox,
                       "Tuple defining the overlap between cell grids for large computations.\n" )
        .def_readonly("cutoff", &t_GenPot::cutoff, "Cutoff of the potential.\n");
    }


  } // namespace python
} // namespace LaDa
