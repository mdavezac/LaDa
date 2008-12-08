//
//  Version: $Id: checkpoints.h 849 2008-11-12 04:45:34Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_DUMMY_H_
#define _LADA_GA_CHECKPOINTS_DUMMY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      namespace Factory
      {
        //! A dummy factory, does nothing.
        class Dummy
        {
          public:
            //! Constructor.
            Dummy() {}
            //! Copy Constructor.
            Dummy( const Dummy& ) {}
            //! dummy.
            template< class T_ARG1 > void operator()( T_ARG1& ) const {}
            //! dummy.
            template< class T_ARG1, class T_ARG2 > void operator()( T_ARG1&, T_ARG2& ) const{}
        };
      }
    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif
