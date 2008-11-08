//
//  Version: $Id: bandgap_stubs.h 816 2008-10-17 01:29:20Z davezac $
//
#ifndef _DARWIN_ELECTRIC_DIPOLE_H_
#define _DARWIN_ELECTRIC_DIPOLE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>
#include <eo/utils/eoRNG.h>

#include <pescan_interface/va.h>
#include <pescan_interface/dipole_elements.h>
#include <crystal/structure.h>
#include <opt/types.h>
#include <mpi/mpi_object.h>
#include <print/stdout.h>

#include "vff.h"

namespace LaDa
{
  namespace GA
  {
    namespace Keepers
    {
      //! Keeps track of oscillator strength
      struct OscStrength
      {
        friend class boost::serialization::access;
        types::t_real osc_strength; //!< Conduction Band Minimum
      
        //! Constructor
        OscStrength() : osc_strength(0) {}
        //! Copy Constructor
        OscStrength(const OscStrength &_c) : osc_strength(_c.osc_strength) {}
        //! Detructor
        ~OscStrength() {};
      
        //! Loads OscStrength::osc_strength from attributes of \a _node.
        bool Load( const TiXmlElement &_node );
        //! Saves OscStrength::osc_strength as attributes of \a _node.
        bool Save( TiXmlElement &_node ) const;
        private:
          //! Serializes a scalar individual.
          template<class Archive> void serialize(Archive & _ar, const unsigned int _version)
            { _ar & osc_strength; }
      };
      //! Prints oscillator strength,
      inline std::ostream& operator<<( std::ostream &_stream, const OscStrength &_o )
        { return _stream << " Transition Dipole: " << _o.osc_strength; }

    }

    namespace OscStrength
    {
      //! \brief GA::Evaluator stub for the oscillator strength functional.
      class Darwin __MPICODE( : public MPI_COMMDEC )
      {
        public:
          //! Constructor and Initializer
          Darwin   ( Crystal::Structure &_s )
                 : __MPICODE( MPI_COMMCOPY( *::LaDa::mpi::main ) __COMMA__ )
                   degeneracy(types::tolerance)  {}
          //! Copy Constructor
          Darwin   ( const Darwin &_b ) 
                 : __MPICODE( MPI_COMMCOPY( _b ) __COMMA__ )
                   degeneracy(_b.degeneracy) {}
          //! Destructor
          ~Darwin() {};

          //! Load Pescan::BandGap from XML
          bool Load( const TiXmlElement &_node );
          //! Computed the oscillator strength, given a bandgap functional.
          void operator()( const Pescan::BandGap &_bg, 
                           const Crystal::Structure& _str, 
                           Keepers::OscStrength& _k )
          {
#           ifndef _NOLAUNCH
              _k.osc_strength = ::LaDa::Pescan::oscillator_strength( _bg, _str, 
                                                               degeneracy, false );
#           else
              _k.osc_strength = rng.uniform();
              Print::out << "osc strength " <<  _k.osc_strength << "\n";
#           endif
          }
#    ifdef _MPI
          //! Sets communicator and suffix for mpi stuff.
          void set_mpi( boost::mpi::communicator *_comm, const std::string &_suffix )
            { MPI_COMMDEC::set_mpi( _comm ); }
#    endif
          //! initializes nothing.
          bool init() { return true; }

        protected:
          //! Allowable Degeneracy.
          types::t_real degeneracy;

      };

    } // namespace ElectricDipoles
  }
} // namespace LaDa

#endif // _PESCAN_H:_
