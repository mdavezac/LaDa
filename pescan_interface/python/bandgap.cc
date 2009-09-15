//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <opt/initial_path.h>

#include "../bandgap.h"
#include "../dipole_elements.h"
#include "bandgap.hpp"

namespace LaDa
{
  namespace Python
  {
    namespace XML
    {
      void Bandgap_from_XML(LaDa::Pescan::BandGap &_type, const std::string &_filename )
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
                       "Could not load Bandgap functional from " + _filename + ".\n" )
      }
    }

    namespace details
    {

      struct pickle_bands : boost::python::pickle_suite
      {
        static boost::python::tuple getinitargs( LaDa::Pescan::Bands const& _w) 
          { return boost::python::make_tuple( _w.vbm, _w.cbm ); }
      };

      std::string print_bands( const LaDa::Pescan::Bands & _b )
      {
        std::ostringstream sstr;
        sstr << "Gap = " << _b.gap() << " = " << _b.cbm << " - " << _b.vbm;
        return sstr.str();
      }

    } // namespace details

    void expose_bands()
    {
      namespace bp = boost::python;
      typedef LaDa::Pescan::Bands t_Bands;
      bp::class_< t_Bands >( "Bands", "Holds vbm and cbm. Used by LaDa.BandGap." )
        .def( bp::init< const LaDa::types::t_real, const LaDa::types::t_real >() )
        .def( bp::init< const t_Bands& >() )
        .def_readwrite( "vbm", &t_Bands::vbm, "Valence Band Maximum." )
        .def_readwrite( "cbm", &t_Bands::cbm, "Conduction Band Minimum." )
        .def( "gap", &t_Bands::gap, "Returns the gap (LaDa.Bands.cbm - LaDa.Bands.vbm)." )
        .def( "__str__", &details::print_bands, "Prints out a Bands object to a string." )
        .def_pickle( details::pickle_bands() );
    }

    void expose_bandgap()
    {
      namespace bp = boost::python;
      typedef LaDa::Pescan::BandGap t_BandGap;
      bp::class_< t_BandGap, bp::bases<LaDa::Pescan::Interface> >
                ( "BandGap", "Computes BandGap. Inherits from LaDa.Escan." )
        .def_readonly( "bands", &t_BandGap :: bands, "BandGap result after calculation." )
        .def_readonly( "vbm_eigs", &t_BandGap :: vbm_eigs, "Eigenvalues from vbm calculation." )
        .def_readonly( "cbm_eigs", &t_BandGap :: cbm_eigs, "Eigenvalues from cbm calculation." )
        .def_readwrite( "eref", &t_BandGap :: Eref,
                        "Reference energies for VBM and CBM calculations.\n"
                        "Is input for Folded-Spectrum calculations, "
                        "and output when performing Full-Diagonalization calculations." )
        .def( "fromXML",  &XML::Bandgap_from_XML, bp::arg("file"),
              "Loads bandgap parameters from an XML file." )
        .def( "evaluate", &t_BandGap::operator(), "Performs a calculation." );
    }

    types::t_real oscillator_strength(Pescan::BandGap const &_bg, types::t_real _d, bool _v)
      { return Pescan::oscillator_strength(_bg, _d, _v); }

    void expose_oscillator_strength()
    {
      namespace bp = boost::python;
      bp::def
      (
        "oscillator_strength",
        &oscillator_strength,
        (
          bp::arg("bandgap_functional"),
          bp::arg("degeneracy") = 0.001,
          bp::arg("verbose") = false
        ),
        "Returns squared norm of the oscillator strength: |<r>|^2."
      );
    }

  } // namespace Python
} // namespace LaDa
