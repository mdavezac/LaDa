//
//  Version: $Id: checkpoints.h 849 2008-11-12 04:45:34Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_APPLY2STORED_H_
#define _LADA_GA_CHECKPOINTS_APPLY2STORED_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      //! \brief Applies a functor to best stored individual.
      //! \details Generally (depending on the overloading of
      //!          Store::Manip::apply_best), this should mean applying to whatever
      //!          optimum exists in Apply2Best::darwin_.store. This functor is not
      //!          permitted to change anything in the stored invididuals.
      template< class T_DARWIN, class T_FUNCTOR >
        class Apply2Stored
        {
          public:
            //! Constructor and Initializer.
            Apply2Stored   ( const T_DARWIN &_darwin, const T_FUNCTOR& _functor )
                         : darwin_(_darwin ), functor_( _functor ) {}
            //! Copy Constructor.
            Apply2Stored   ( const Apply2Stored &_c )
                         : darwin_(_c.darwin_ ), functor_( _c.functor_ ) {}
       
            //! Functor. Reroutes calls to Store::Manip::apply_best().
            void operator()( bool ) { darwin_.store->apply_all( functor ); }

          protected:
            //! The store functor.
            const T_DARWIN &darwin_;
            //! The application functor.
            const T_FUNCTOR &functor_;
        };
      
      //! Returns an Apply2Best functor.
      template< class T_STORE, class T_FUNCTOR >
        Apply2Stored<T_DARWIN, T_FUNCTOR> apply2best( const T_DARWIN& _darwin, const T_FUNCTOR& _functor )
          { return Apply2Stored<T_DARWIN, T_FUNCTOR> Apply2Stored( _darwin, _functor ); }
  
    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif
