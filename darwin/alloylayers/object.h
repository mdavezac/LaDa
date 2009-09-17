//
//  Version: $Id$
//
#ifndef  _DARWIN_ALLOYLAYERS_OBJECT_H_
#define  _DARWIN_ALLOYLAYERS_OBJECT_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <string>

#include<opt/types.h>
#include<vff/layered.h>

#include "../print_callbacks.h"
#include "../bitstring/object.h"
#include "../static_translate.h"
#include "../effective_mass.h"
#include "../bandgap_stubs.h"
#include "../electric_dipole.h"

#include "concentration.h"

namespace LaDa
{
  //! \cond
  namespace Crystal
  {
    class Structure;
  }
  //! \endcond

  namespace GA
  {
    namespace AlloyLayers
    {
      struct Object : public LaDa::GA::BitString::Object<>, 
                      public LaDa::GA::Keepers::BandGap,
                      public LaDa::GA::Keepers::OscStrength,
                      public LaDa::GA::Keepers::Concentration,
                      public LaDa::GA::Keepers::eMass,
                      public LaDa::GA::Keepers::hMass,
                      public LaDa::GA::PrintSignal< Object >
      {
        friend class boost::serialization::access;
        //! The type of the BitString container
        typedef LaDa::GA::BitString::Object<> :: t_Container t_Container;
        //! Constructor
        Object() : LaDa::GA::BitString::Object<>(),
                   LaDa::GA::Keepers::BandGap(),
                   LaDa::GA::Keepers::OscStrength(),
                   LaDa::GA::Keepers::Concentration(),
                   LaDa::GA::Keepers::eMass(),
                   LaDa::GA::Keepers::hMass() {}
        //! Copy Constructor
        Object   (const Object &_c)
               : LaDa::GA::BitString::Object<>(_c),
                 LaDa::GA::Keepers::BandGap(_c),
                 LaDa::GA::Keepers::OscStrength(_c),
                 LaDa::GA::Keepers::Concentration( _c ),
                 LaDa::GA::Keepers::eMass( _c ),
                 LaDa::GA::Keepers::hMass( _c ) {};
        //! Loads from \a _node.
        bool Load( const TiXmlElement &_node )
          { return     LaDa::GA::Keepers::BandGap::Load( _node )
                   and LaDa::GA::Keepers::OscStrength::Load( _node )
                   and LaDa::GA::Keepers::Concentration::Load( _node )
                   and LaDa::GA::Keepers::eMass::Load( _node )
                   and LaDa::GA::Keepers::hMass::Load( _node ); }
        //! Saves to \a _node.
        bool Save( TiXmlElement &_node ) const
          { return     LaDa::GA::Keepers::BandGap::Save( _node )
                   and LaDa::GA::Keepers::OscStrength::Save( _node ) 
                   and LaDa::GA::Keepers::Concentration::Save( _node )
                   and LaDa::GA::Keepers::eMass::Save( _node )
                   and LaDa::GA::Keepers::hMass::Save( _node ); }
        //! Destructor
        virtual ~Object() {};
        //! \brief Comparison operator.
        //! \details Checks symmetric equivalents.
        bool operator==( Object const &_c ) const;
        //! Chooses random symmetry from all possible symmetries.
        bool random_symmetry();
        private:
          //! Serializes a scalar individual.
          template<class Archive>
            void serialize(Archive & _ar, const unsigned int _version)
            {
              _ar & boost::serialization::base_object< LaDa::GA::BitString::Object<> >( *this ); 
              _ar & boost::serialization::base_object< LaDa::GA::Keepers::BandGap >( *this ); 
              _ar & boost::serialization::base_object< LaDa::GA::Keepers::OscStrength >( *this ); 
              _ar & boost::serialization::base_object< LaDa::GA::Keepers::Concentration >( *this ); 
              _ar & boost::serialization::base_object< LaDa::GA::Keepers::eMass >( *this ); 
              _ar & boost::serialization::base_object< LaDa::GA::Keepers::hMass >( *this ); 
            }
      };

    }
  }
} // namespace LaDa
#endif
