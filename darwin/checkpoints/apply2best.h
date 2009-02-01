//
//  Version: $Id: checkpoints.h 849 2008-11-12 04:45:34Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_APPLY2BEST_H_
#define _LADA_GA_CHECKPOINTS_APPLY2BEST_H_

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
        class Apply2Best 
        {
          public:
            //! Constructor and Initializer.
            Apply2Best   ( const T_DARWIN &_darwin, const T_FUNCTOR& _functor )
                       : darwin_(_darwin ), functor_( _functor ) {}
            //! Copy Constructor.
            Apply2Best   ( const Apply2Best &_c )
                       : darwin_(_c.darwin_ ), functor_( _c.functor_ ) {}
       
            //! Functor. Reroutes calls to Store::Manip::apply_best().
            bool operator()( bool )
             { darwin_.store->apply_best( functor_ ); return true; }

          protected:
            //! The darwin object.
            const T_DARWIN &darwin_;
            //! The application functor.
            const T_FUNCTOR &functor_;
        };
      
      //! Returns an Apply2Best functor.
      template< class T_DARWIN, class T_FUNCTOR >
        Apply2Best<T_DARWIN, T_FUNCTOR> apply2best( const T_DARWIN& _darwin, const T_FUNCTOR& _functor )
          { return Apply2Best<T_DARWIN, T_FUNCTOR>( _darwin, _functor ); }
  
    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif
