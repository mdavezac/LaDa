//
//  Version: $Id$
//
#ifndef  _DARWIN_ALLOYLAYERS_OBJECT_H_
#define  _DARWIN_ALLOYLAYERS_OBJECT_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

namespace GA
{
  namespace AlloyLayers
  {
    struct Object : public Bitstring::Object<>, 
                    public ::GA::Keepers::BandGap,
                    public ::GA::Keepers::OscStrength
    {
      //! The type of the BitString container
      typedef ::BitString::Object<> :: t_Container t_Container;
      //! Constructor
      Object() : ::BitString::Object<>(), ::GA::Keepers::BandGap(),
                 ::GA::Keepers::OscStrength() {}
      //! Copy Constructor
      Object   (const Object &_c)
             : ::BitString::Object<>(_c), ::GA::Keepers::BandGap(_c),
               ::GA::Keepers::OscStrength(_c) {};
      //! Loads from \a _node.
      bool Load( const TiXmlElement &_node )
        { return     ::GA::Keepers::BandGap::Load( _node ) 
                 and ::GA::Keepers::OscStrength::Load( _node ); }
      //! Saves to \a _node.
      bool Save( TiXmlElement &_node ) const
        { return     ::GA::Keepers::BandGap::Save( _node ) 
                 and ::GA::Keepers::OscStrength::Save( _node ); }
      //! Destructor
      virtual ~Object() {};
      //! Serializes a scalar individual.
      template<class Archive>
        void serialize(Archive & _ar, const unsigned int _version)
        {
          _ar & boost::serialization::base_object< ::BitString::Object<> >( *this ); 
          _ar & boost::serialization::base_object< ::GA::Keepers::BandGap >( *this ); 
          _ar & boost::serialization::base_object<::GA::Keepers::OscStrength>( *this ); 
        }
    };

  }
}

#endif
