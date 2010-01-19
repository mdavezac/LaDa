//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/def.hpp>
#include <boost/python/str.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/data_members.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <pyublas/numpy.hpp>

#include <opt/initial_path.h>

#include "../interface.h"
#include "escan.hpp"


namespace LaDa
{
  namespace Python
  {
    namespace XML
    {
      template<class T_TYPE>
        void Escan_from_XML(T_TYPE &_type, const std::string &_filename )
        {
          TiXmlDocument doc( _filename ); 
          TiXmlHandle docHandle( &doc ); 
          if( not opt::InitialPath::is_initialized() )
            opt::InitialPath::init();
        
          __DOASSERT( not doc.LoadFile(), 
                         doc.ErrorDesc() << "\n"  
                      << "Could not load input file " << _filename  
                      << ".\nAborting.\n" ) 
          __DOASSERT( not docHandle.FirstChild("Job").Element(),
                      "Could not find <Job> tag in " << _filename << ".\n" )
         
          __DOASSERT( not _type.Load( *docHandle.FirstChild("Job").Element() ),
                         "Could not load Pescan functional from " + _filename + ".\n" )
        }
    }

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
#   undef LADA_PARAMS

    typedef LaDa::Pescan::Interface::GenPot::t_MeshTuple t_MeshTuple;
    
    boost::python::tuple get_impl( t_MeshTuple const& _mesh )
    {
      namespace bt = boost::tuples;
      namespace bp = boost::python;
      return bp::make_tuple( bt::get<0>(_mesh), bt::get<1>(_mesh), bt::get<2>(_mesh) );
    }
    void set_impl( boost::python::tuple const &_tuple, t_MeshTuple & _mesh )
    {
      namespace bt = boost::tuples;
      namespace bp = boost::python;
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
    
    boost::python::tuple get_mesh( LaDa::Pescan::Interface::GenPot const &_genpot )
      { return get_impl(_genpot.mesh); }
    void set_mesh( LaDa::Pescan::Interface::GenPot &_genpot, boost::python::tuple const &_t )
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
    boost::python::tuple get_mcell( LaDa::Pescan::Interface::GenPot const &_genpot )
      { return get_impl(_genpot.multiple_cell); }
    void set_mcell( LaDa::Pescan::Interface::GenPot &_genpot, boost::python::tuple const &_t )
      { set_impl(_t, _genpot.multiple_cell); }
    boost::python::tuple get_sbox( LaDa::Pescan::Interface::GenPot const &_genpot )
      { return get_impl(_genpot.small_box); }
    void set_sbox( LaDa::Pescan::Interface::GenPot &_genpot, boost::python::tuple const &_t )
      { set_impl(_t, _genpot.small_box); }

    boost::python::str get_dir( LaDa::Pescan::Interface const &_self )
      { return boost::python::str( _self.get_dirname().string() ); }
    void set_dir( LaDa::Pescan::Interface &_self, boost::python::str const &_str )
    {
      std::string const string = boost::python::extract<std::string>(_str);
      return _self.set_dirname(string); 
    }

    pyublas::numpy_vector<types::t_real> get_eigenvalues(Pescan::Interface const& _interface)
    {
      pyublas::numpy_vector<types::t_real> a(_interface.eigenvalues.size());
      std::copy(_interface.eigenvalues.begin(), _interface.eigenvalues.end(), a.begin());
      return a;
    }

    void expose_escan()
    {
      typedef LaDa::Pescan::Interface t_Escan;
      namespace bp = boost::python;
      bp::enum_<t_Escan::t_method>( "method", "Diagonalisation method." )
        .value( "folded", t_Escan::FOLDED_SPECTRUM )
        .value( "full_diagonalization", t_Escan::ALL_ELECTRON )
        .export_values();
      bp::enum_<t_Escan::Escan::t_potential>( "potential", "Type of the hamiltonian" )
        .value( "local", t_Escan::Escan::LOCAL )
        .value( "nonlocal", t_Escan::Escan::NONLOCAL )
        .value( "spinorbit", t_Escan::Escan::SPINORBIT )
        .export_values();

      bp::class_< t_Escan >( "Escan", "Wrapper around the nanopse-escan functional." ) 
        .def( bp::init< t_Escan& >() )
        .add_property
        (
           "directory", &get_dir, &set_dir, 
           "Directory where calculations are performed.\n"
        )
        .add_property
        ( 
          "scale", &t_Escan::get_scale, &t_Escan::set_scale,
          "Internal escan scale.\n\n"
          "Prior to calculation, set as C{escan.scale = structure} where \"structure\" is a "
          "L{sStructure<lada.crystal.sStructure>} or L{Structure<lada.crystal.Structure>} object.\n"
        )
#       define LADA_PARAMS(name, docstring)\
          .add_property(#name, &get ## name, &set ## name, docstring)
          LADA_PARAMS(_input_filename, "Input filename for pescan.")
          LADA_PARAMS(_output_filename, "Output filename for pescan.")
          LADA_PARAMS(_rcut, "Real-space cutoff.")
          LADA_PARAMS(_smoothness, "Smoothness factor of plane-wave cutoff.")
          LADA_PARAMS(_kinscal, "Kinetic scaling parameter.")
#       undef LADA_PARAMS
#       define LADA_PARAMS(name, docstring) \
          .add_property(#name, &get_ ## name, &set_ ## name, docstring)
          LADA_PARAMS(vff_inputfile, "Structure input file.")
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
           bp::make_function(&set_kpoint, bp::return_value_policy<bp::return_by_value>()),
           "K-point, in units of M{2S{pi}/a}, at which to perform calculation."
        )
        .def_readwrite( "destroy_directory", &t_Escan::do_destroy_dir,
                        "If true, directory where calculations are carried out is destroyed "
                        "at the end of a calculation." )
        .add_property( "eigenvalues", &get_eigenvalues, 
                       "Computed eigenvalues as a read-only numpy array of real values." )
        .def_readwrite("genpot", &t_Escan::genpot, "Parameters of the atomic potentials.")
        .def_readwrite("verbose", &t_Escan::verbose, "Verbose pescan output on true.")
        .def( "fromXML",  &XML::Escan_from_XML<t_Escan>, bp::arg("file"),
              "Loads escan parameters from an XML file." )
        .def( "__call__", &t_Escan::operator(), "Performs a calculation." )
        .def( "set_mpi", &t_Escan::set_mpi, "Sets the boost.mpi communicator." );

      bp::def( "nb_valence_states", &LaDa::Pescan::nb_valence_states, bp::arg("structure"),
               "Returns the number of valence states in a structure." );
    }


    void expose_genpot()
    {

      typedef LaDa::Pescan::Interface::GenPot t_GenPot;
      namespace bp = boost::python;
      bp::class_< t_GenPot >( "GenPot", "Wrapper around genpot parameters.", bp::no_init ) 
        .add_property( "mesh", &get_mesh, &set_mesh ) 
        .add_property( "multiple_cell", &get_mcell, &set_mcell ) 
        .add_property( "small_box", &get_sbox, &set_sbox )
        .def_readonly("cutoff", &t_GenPot::cutoff);
    }


//   void expose_escan_parameters()
//   {
//     bp::class_< t_Parameters >( "EscanParameters", bp::init< t_Parameters&>() ) 
//       .add_property( EXPOSE_FILENAME( input_filename ), 
//                      "Filename where to write the input of nanopse-escan input." 
//                      "You should not need to touch this." )
//       .add_property( EXPOSE_FILENAME( output_filename ),
//                      "Filename of the output of nanopse-escan."
//                      "You should not need to touch this." )
//       .def_readwrite( "Eref", &t_Parameters::Eref,
//                       "Reference energy for folded-spectrum calculation." )
//       .def_readwrite( "nbstates", &t_Parameters::nbstates,
//                       "Number of states to compute." )
//       .def_readwrite( "kpoint", &t_Parameters::kpoint,
//                       "kpoint in units of 2pi/a, where a is the unit-cell parameter."
//                       "If a cell is not cubic, then check nanopse-escan documentation." )
//       .def_readwrite( "method", &t_Parameters::method,
//                       "Calculation method: \"folded\" or \"full_diagonalization\"." )
//       .def_readonly( "potential", &t_Parameters::potential,
//                      "Type of potential: \"local\", \"nonlocal\", or \"spinorbit\"." );
//   }
  } // namespace Python
} // namespace LaDa
