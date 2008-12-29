//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <pescan_interface/interface.h>
#include <opt/initial_path.h>

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
                         "Could not load Vff functional from " + _filename + ".\n" )
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

      template< class T_TYPE, size_t LINESIZE >
      std::string print_vector( const std::vector<T_TYPE> & _container )
      {
        std::ostringstream sstr;
        size_t i(0);
        foreach( const T_TYPE& value, _container )
        {
          if( i % LINESIZE == 0 and i != 0 ) sstr << "\n";
          sstr << value << " ";
          ++i;
        }
        return sstr.str();
      }

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
        .add_property( EXPOSE_FILENAME( input_filename ) )
        .add_property( EXPOSE_FILENAME( output_filename ) )
        .def_readwrite( "Eref", &t_Parameters::Eref )
        .def_readwrite( "nbstates", &t_Parameters::nbstates )
        .def_readwrite( "kpoint", &t_Parameters::kpoint )
        .def_readwrite( "method", &t_Parameters::method )
        .def_readonly("potential", &t_Parameters::potential );
    }

    void expose_escan()
    {
      typedef LaDa::Pescan::Interface t_Escan;
      namespace bp = boost::python;
      bp::class_< t_Escan >( "Escan" ) 
        .def( bp::init< t_Escan& >() )
        .add_property( EXPOSE_FILENAME( vff_inputfile ) )
        .add_property( EXPOSE_FILENAME( directory ) )
        .add_property( "scale", &t_Escan::get_scale, &t_Escan::set_scale )
        .def_readwrite( "parameters", &t_Escan::escan )
        .def_readwrite( "destroy_directory", &t_Escan::do_destroy_dir )
        .def_readonly( "eigenvalues", &t_Escan::eigenvalues )
        .def( "fromXML",  &XML::Escan_from_XML<t_Escan> )
        .def( "run", &t_Escan::operator() )
        .def( "set_mpi", &t_Escan::set_mpi );
      // Register eigenvalues type.
      bp::class_< std::vector<types::t_real> >("VectorReals")
       .def(bp::vector_indexing_suite< std::vector<types::t_real> >())
       .def("__str__", &details::print_vector< types::t_real, 5 > );

      bp::def( "nb_valence_states", &LaDa::Pescan::nb_valence_states );
    }

#   undef EXPOSE_FILENAME
  } // namespace Python
} // namespace LaDa
