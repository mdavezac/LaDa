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

#include "../bitstring.h"
#include "policies.h"

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
      struct Object : public LaDa::BitString::Object<>, 
                      public LaDa::GA::Keepers::BandGap,
                   //  public LaDa::GA::Keepers::OscStrength,
                      public LaDa::GA::AlloyLayers::PrintSignal< Object >
      {
        friend class boost::serialization::access;
        //! The type of the BitString container
        typedef LaDa::BitString::Object<> :: t_Container t_Container;
        //! Constructor
        Object() : LaDa::BitString::Object<>(), LaDa::GA::Keepers::BandGap() {}
//                  LaDa::GA::Keepers::OscStrength() {}
        //! Copy Constructor
        Object   (const Object &_c)
               : LaDa::BitString::Object<>(_c), LaDa::GA::Keepers::BandGap(_c) {}
//                LaDa::GA::Keepers::OscStrength(_c) {};
        //! Loads from \a _node.
        bool Load( const TiXmlElement &_node )
          { return     LaDa::GA::Keepers::BandGap::Load( _node ); } 
//                  and LaDa::GA::Keepers::OscStrength::Load( _node ); }
        //! Saves to \a _node.
        bool Save( TiXmlElement &_node ) const
          { return     LaDa::GA::Keepers::BandGap::Save( _node ); }
//                  and LaDa::GA::Keepers::OscStrength::Save( _node ); }
        //! Destructor
        virtual ~Object() {};
        private:
          //! Serializes a scalar individual.
          template<class Archive>
            void serialize(Archive & _ar, const unsigned int _version)
            {
              _ar & boost::serialization::base_object< LaDa::BitString::Object<> >( *this ); 
              _ar & boost::serialization::base_object< LaDa::GA::Keepers::BandGap >( *this ); 
//             _ar & boost::serialization::base_object< LaDa::GA::Keepers::OscStrength >( *this ); 
            }
      };

      //! \brief Old-style translation.
      //! \todo remove this function and replace it with translate.
      inline void operator<<( std::string &_str, const Object& _o )
        { Translate< Object > :: translate( _o, _str ); }
      //! \brief Old-style translation.
      //! \todo remove this function and replace it with translate.
      inline void operator<<( Object& _o, const std::string &_str )
        { Translate< Object > :: translate( _str, _o ); }
      //! \brief Old-style translation.
      //! \todo remove this function and replace it with translate.
      inline void operator<<( Crystal::Structure &_str, const Object& _o )
        { Translate< Object > :: translate( _o, _str ); }
      //! \brief Old-style translation.
      //! \todo remove this function and replace it with translate.
      inline void operator<<( Object& _o, const Crystal::Structure &_str )
        { Translate< Object > :: translate( _str, _o ); }
    }
  }
} // namespace LaDa
#endif
