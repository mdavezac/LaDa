//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/def.hpp>
#include <boost/python/enum.hpp>

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

    namespace details
    {

#     ifdef FILETOPATH_GETSETTERS
#       error Macro FILETOPATH_GETSETTERS already exists.
#     endif
#     define FILETOPATH_GETSETTERS( _C_, _A_, _B_ ) \
        const std::string get_ ## _B_ ( const _C_ &_a )  { return _a._A_.string(); } \
        void set_ ## _B_ ( _C_ &_a, const std::string &_b )  { _a._A_ = _b; }
#     ifdef EXPOSE_FILENAME
#       error Macro EXPOSE_FILENAME already exists.
#     endif
      
      FILETOPATH_GETSETTERS( LaDa::Pescan::Interface::Escan, filename, input_filename )
      FILETOPATH_GETSETTERS( LaDa::Pescan::Interface::Escan, output, output_filename )
      FILETOPATH_GETSETTERS( LaDa::Pescan::Interface, atom_input, vff_inputfile )
      const std::string get_directory ( const LaDa::Pescan::Interface &_a ) 
        { return _a.get_dirname().string(); } 
      void set_directory ( LaDa::Pescan::Interface &_a, const std::string &_b ) 
        { _a.set_dirname( _b ); }

    } // namespace details

#   undef FILETOPATH_GETSETTERS
#   define EXPOSE_FILENAME( _B_ ) \
     #_B_,  &details::get_ ## _B_, &details::set_ ## _B_

    void expose_escan_parameters()
    {
      typedef LaDa::Pescan::Interface t_Escan;
      typedef LaDa::Pescan::Interface::Escan t_Parameters;
      namespace bp = boost::python;
      bp::enum_< t_Escan::t_method >( "method" )
        .value( "folded", t_Escan::FOLDED_SPECTRUM )
        .value( "full_diagonalization", t_Escan::ALL_ELECTRON )
        .export_values();
      bp::enum_< t_Parameters::t_potential >( "potential" )
        .value( "local", t_Parameters::LOCAL )
        .value( "nonlocal", t_Parameters::NONLOCAL )
        .value( "spinorbit", t_Parameters::SPINORBIT )
        .export_values();

      bp::class_< t_Parameters >( "EscanParameters", bp::init< t_Parameters&>() ) 
        .add_property( EXPOSE_FILENAME( input_filename ), 
                       "Filename where to write the input of nanopse-escan input." 
                       "You should not need to touch this." )
        .add_property( EXPOSE_FILENAME( output_filename ),
                       "Filename of the output of nanopse-escan."
                       "You should not need to touch this." )
        .def_readwrite( "Eref", &t_Parameters::Eref,
                        "Reference energy for folded-spectrum calculation." )
        .def_readwrite( "nbstates", &t_Parameters::nbstates,
                        "Number of states to compute." )
        .def_readwrite( "kpoint", &t_Parameters::kpoint,
                        "kpoint in units of 2pi/a, where a is the unit-cell parameter."
                        "If a cell is not cubic, then check nanopse-escan documentation." )
        .def_readwrite( "method", &t_Parameters::method,
                        "Calculation method: \"folded\" or \"full_diagonalization\"." )
        .def_readonly( "potential", &t_Parameters::potential,
                       "Type of potential: \"local\", \"nonlocal\", or \"spinorbit\"." );
    }

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

    void expose_escan()
    {
      typedef LaDa::Pescan::Interface t_Escan;
      namespace bp = boost::python;
      bp::class_< t_Escan >( "Escan", "Wrapper around the nanopse-escan functional." ) 
        .def( bp::init< t_Escan& >() )
        .add_property( EXPOSE_FILENAME( vff_inputfile ), "Structure input file." )
        .add_property( EXPOSE_FILENAME( directory ), "Directory where to perform calculations." )
        .add_property( "scale", &t_Escan::get_scale, &t_Escan::set_scale,
                       "Internal escan scale. Prior to calculation, "
                       "set as \"escan.scale = structure\" "
                       "where \"structure\" is a LaDa.Structure object." )
        .def_readwrite( "parameters", &t_Escan::escan,
                        "EscanParameters object." )
        .def_readwrite( "destroy_directory", &t_Escan::do_destroy_dir,
                        "If true, directory where calculations are carried out is destroyed "
                        "at the end of a calculation." )
        .def_readonly( "eigenvalues", &t_Escan::eigenvalues,
                       "Computed eigenvalues." )
        .def_readwrite("genpot", &t_Escan::genpot)
        .def( "fromXML",  &XML::Escan_from_XML<t_Escan>, bp::arg("file"),
              "Loads escan parameters from an XML file." )
        .def( "run", &t_Escan::operator(), "Performs a calculation." )
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

#   undef EXPOSE_FILENAME
  } // namespace Python
} // namespace LaDa
