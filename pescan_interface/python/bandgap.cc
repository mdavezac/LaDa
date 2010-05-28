#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <opt/initial_path.h>

#include "../bandgap.h"
#include "../dipole_elements.h"
#include "bandgap.hpp"
#include "escan.hpp"
#include <vff/python/vff.hpp>

namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    void Bandgap_from_XML(LaDa::Pescan::BandGap &_type, const std::string &_filename )
    {
      TiXmlDocument doc( _filename ); 
      TiXmlHandle docHandle( &doc ); 
      if( not opt::InitialPath::is_initialized() ) opt::InitialPath::init();
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
                 "Encountered error while loading bandgap functional from xml file "
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
      Bandgap_from_XML(*result, _filename); 
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
        Bandgap_from_XML(*result, _filename); 
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
    

    namespace details
    {

      struct pickle_bands : bp::pickle_suite
      {
        static bp::tuple getinitargs( LaDa::Pescan::Bands const& _w) 
          { return bp::make_tuple( _w.vbm, _w.cbm ); }
      };

      std::string print_bands( const LaDa::Pescan::Bands & _b )
      {
        std::ostringstream sstr;
        sstr << "Gap = " << _b.gap() << " = " << _b.cbm << " - " << _b.vbm;
        return sstr.str();
      }


    } // namespace details

    template<class T>
      Pescan::Bands __call2__(Pescan::BandGap &_interface,
                              T &_vff, Crystal::TStructure<std::string> const &_str)
      {
        // Checks that directory exists.
        namespace bfs = boost::filesystem;
        Pescan::Interface::t_Path rootdir = opt::InitialPath::path()/_interface.get_dirname();
        if( _interface.comm().rank() == 0 and (not python::create_directory(rootdir)) )
        {
          PyErr_SetString( PyExc_IOError,  
                           ("Could not find/create directory: " + rootdir.string()).c_str() );
          bp::throw_error_already_set();
          return _interface.bands;
        }
        _interface.comm().barrier();
        // creates atom input
        Pescan::Interface::t_Path const path = _interface.atom_input;
        Pescan::Interface::t_Path filename = path.filename();
        if( filename.empty() ) filename = "atom_input";
        boost::mpi::communicator world;
        std::ostringstream sstr; sstr << world.rank();
        _interface.atom_input = rootdir / filename.replace_extension(sstr.str());
        print_escan_input(_vff, _interface.atom_input.string(), _str);
        // finally, computes.
        try { _interface(_str); }
        catch(...)
        {
          PyErr_SetString(PyExc_RuntimeError, "Something went wrong in escan.\n");
          bp::throw_error_already_set();
          return _interface.bands;
        }
        _interface.atom_input = path;
        return _interface.bands;
      }

    void expose_bands()
    {
      typedef LaDa::Pescan::Bands t_Bands;
      bp::class_< t_Bands >( "Bands", "Holds vbm and cbm. Used by L{BandGap}.\n\n" )
        .def( bp::init< const LaDa::types::t_real, const LaDa::types::t_real >() )
        .def( bp::init< const t_Bands& >() )
        .def_readwrite( "vbm", &t_Bands::vbm, "Valence Band Maximum." )
        .def_readwrite( "cbm", &t_Bands::cbm, "Conduction Band Minimum." )
        .add_property( "gap", &t_Bands::gap, "Returns the gap (LaDa.Bands.cbm - LaDa.Bands.vbm)." )
        .def( "__str__", &details::print_bands, "Prints out a Bands object to a string." )
        .def_pickle( details::pickle_bands() );
    }

    void expose_bandgap()
    {
      typedef LaDa::Pescan::BandGap t_BandGap;
      bp::class_< t_BandGap, bp::bases<LaDa::Pescan::Interface> >
      (
        "BandGap",
        "Computes BandGap.\n\n"
        "  - No argument.\n"\
        "  - A single string argument representing the path to an XML input file.\n" 
        "  - A single boost.mpi communicator.\n" 
        "  - A string argument (see above), followed by a boost.mpi communicator.\n\n"  
        "If compiled without mpi, including the communicator will have no effect.\n"
      )
        .def( bp::init<t_BandGap const&>() )
        .def( "__init__", bp::make_constructor(&cload<t_BandGap>) )\
        .def( "__init__", bp::make_constructor(&mpi_create<t_BandGap>) )\
        .def( "__init__", bp::make_constructor(&mpi_create2<t_BandGap>) )\
        .def_readonly( "bands", &t_BandGap :: bands, "BandGap result after calculation." )
        .def_readonly( "vbm_eigs", &t_BandGap :: vbm_eigs, "Eigenvalues from vbm calculation." )
        .def_readonly( "cbm_eigs", &t_BandGap :: cbm_eigs, "Eigenvalues from cbm calculation." )
        .def_readwrite( "eref", &t_BandGap :: Eref,
                        "Reference energies for VBM and CBM calculations.\n"
                        "Is input for Folded-Spectrum calculations, "
                        "and output when performing Full-Diagonalization calculations." )
        .def( "fromXML",  &Bandgap_from_XML, bp::arg("file"),
              "Loads bandgap parameters from an XML file." )
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
          "Returns a numpy vector also available as self.L{eigenvalues}. As "
          "such, it will change from one calculation to the other. Copy if you "
          "want to store the results.\n" 
          "@param vff: A vff functional.\n"
          "@type vff: L{vff.Vff} or L{vff.LayeredVff}\n"
          "@param structure: The structure for which to perform calculations.\n"
          "@type structure: L{crystal.Structure}\n",
          bp::with_custodian_and_ward_postcall<1,0>() 
        );
    }
  } // namespace Python
} // namespace LaDa
