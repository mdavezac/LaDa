//
//  Version: $Id$
//
#ifndef  _DARWIN_GROUNDSTATES_OBJECT_H_
#define  _DARWIN_GROUNDSTATES_OBJECT_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <string>

#include<opt/types.h>
#include<vff/layered.h>

#include "../bitstring/object.h"
#include "assign.h"

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
    namespace GroundStates
    {
      struct Object : public LaDa::BitString::Object<>, 
                      public LaDa::GA::Keepers::CE
      {
        friend class boost::serialization::access;
        //! The type of the BitString container
        typedef LaDa::BitString::Object<> :: t_Container t_Container;
        //! Constructor
        Object() : LaDa::BitString::Object<>(), LaDa::GA::Keepers::CE() {}
        //! Copy Constructor
        Object   (const Object &_c)
               : LaDa::BitString::Object<>(_c), LaDa::GA::Keepers::CE(_c) {}
        //! Loads from \a _node.
        bool Load( const TiXmlElement &_node )
          { return  LaDa::GA::Keepers::CE::Load( _node ); } 
        //! Saves to \a _node.
        bool Save( TiXmlElement &_node ) const
          { return  LaDa::GA::Keepers::CE::Save( _node ); }
        //! Destructor
        virtual ~Object() {};
        private:
          //! Serializes a scalar individual.
          template<class Archive>
            void serialize(Archive & _ar, const unsigned int _version)
            {
              _ar & boost::serialization::base_object< LaDa::BitString::Object<> >( *this ); 
              _ar & boost::serialization::base_object< LaDa::GA::Keepers::CE >( *this ); 
            }
      };

      //! Dumps object to a string.
      inline std::ostream& operator<<( std::ostream& _stream, const Object& _o )
      {
        foreach( types::t_real var, _o.Container() )
            _stream << ( var > 0 ? "1": "0" );
        return _stream << "  -  Energy: " << _o.energy;
      }

    }
  }
} // namespace LaDa
#endif
