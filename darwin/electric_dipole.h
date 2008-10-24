//
//  Version: $Id: bandgap_stubs.h 816 2008-10-17 01:29:20Z davezac $
//
#ifndef _DARWIN_BANDGAP_STUBS_H_
#define _DARWIN_BANDGAP_STUBS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

#include <pescan_interface/va.h>
#include <crystal/structure.h>
#include <opt/types.h>
#include <mpi/mpi_object.h>

#include "vff.h"

namespace GA
{
  namespace Keepers
  {
    //! Keeps track of oscillator strength
    struct OscStrength
    {
      types::t_real osc_strength; //!< Conduction Band Minimum
    
      //! Constructor
      OscStrength() : osc_strength(0) {}
      //! Copy Constructor
      OscStrength(const OscStrength &_c) : osc_strength(_c.osc_strength) {}
      //! Detructor
      ~OscStrength() {};
    
      //! Loads BandGap::vbm and BandGap::cvm from attributes of \a _node.
      bool Load( const TiXmlElement &_node );
      //! Saves BandGap::vbm and BandGap::cvm as attributes of \a _node.
      bool Save( TiXmlElement &_node ) const;
      //! Serializes a scalar individual.
      template<class Archive> void serialize(Archive & _ar, const unsigned int _version)
        { _ar & osc_strength; }
    };
  }

  namespace OscStrength
  {
    //! \brief GA::Evaluator stub for the oscillator strength functional.
    class Darwin __MPICODE( : public MPI_COMMDEC )
    {
      protected:
        //! \brief Structure on which to perform calculations
        //! \details The structure which is referenced should be updated prior to
        //!          calling Darwin::operator()(). Although much of the
        //!          information should come from Darwin::atomicconfig, such
        //!          things as the number of electrons (for all-electron
        //!          calculations) does not.
        Crystal::Structure &structure;
        //! Allowable Degeneracy.
        types::t_real degeneracy;

      public:
        //! Constructor and Initializer
        Darwin   ( Crystal::Structure &_s )
               : __MPICODE( MPI_COMMCOPY( *::mpi::main ) __COMMA__ )
                 structure(_s), degeneracy(types::t_real) 
                 __MPICODE( __COMMA__ suffix("") ) {}
        //! Copy Constructor
        Darwin   ( const Darwin &_b ) 
               : __MPICODE( MPI_COMMCOPY( _b ) __COMMA__ )
                 structure(_b.structure), degeneracy(_b.degeneracy) 
                 __MPICODE( __COMMA__ suffix( _b.suffix ) ) {}
        //! Destructor
        ~Darwin() {};

        //! Load Pescan::BandGap from XML
        bool Load( const TiXmlElement &_node )
        {
          const TiXmlElement *parent;
          std::string str;
      
          str = _node.Value();
          if ( str.compare("Functional" ) != 0 )
            parent = _node.FirstChildElement("Functional");
          else
            parent = &_node;
      
          while (parent)
          {
            str = "";
            if ( parent->Attribute( "type" )  )
              str = parent->Attribute("type");
            if ( str.compare("oscillator strength" ) == 0 ) break;
            parent = parent->NextSiblingElement("Functional");
          }
          if( not parent ) return false;
          if( not parent->Attribute( "degeneracy" ) ) return false;
          degeneracy = boost::lexical_cast< types::t_real >
                                          ( parent->Attribute("degeneracy") );
          return true;
        }
        //! Computed the oscillator strength, given a bandgap functional.
        void operator()( const ::Pescan::BandGap &_bg, 
                         const Crystal::Structure& _str, 
                         ::GA::Keepers::OscStrength& _k );
        {
          _k.osc_strength = ::Pescan::oscillator_strength( _structure, _bg,
                                                           degeneracy, false );
        }
  #ifdef _MPI
        //! Sets communicator and suffix for mpi stuff.
        void set_mpi( boost::mpi::communicator *_comm, const std::string &_suffix )
          { MPI_COMMDEC::set_mpi( _comm ); }
  #endif
        //! initializes the vff part
        bool init() {}

    };

  } // namespace ElectricDipoles
}

#endif // _PESCAN_H_
