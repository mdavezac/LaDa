//
//  Version: $Id$
//
#ifndef _DARWIN_GROUNDSTATES_ASSIGN_H_
#define _DARWIN_GROUNDSTATES_ASSIGN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa
{
  namespace GA
  {
    namespace GroundStates
    {
      //! Small policy to assign ce result to single-value quantity.
      template< class T_OBJECT, class T_QUANTITIES >
        class AssignCE
        {
          public:
            //! Type of the object.
            typedef T_OBJECT t_Object;
            //! Type of the quantities.
            typedef T_QUANTITIES t_Quantities;
       
            //! Assigns a value from object to a quantity.
            void assign( const t_Object& _o, t_Quantities &_q ) const
              { _q = _o.::LaDa::Ga::Keepers::CE::energy; }
        };

    }
  }
}
#endif
