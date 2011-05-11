//
//  Version: $Id: current.h 1163 2009-06-07 22:15:19Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_CURRENT_H_
#define _LADA_LOADNSAVE_TRANSFORMS_CURRENT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      //! Returns current expression.
      struct current : proto::callable
      {
        //! Resulting type.
        template< class T_SIGNA > struct result;
        //! Specialization for calls.
        template< class T_THIS, class T_EXPR >
          struct result< T_THIS(T_EXPR const&) >
          {
            //! Returns a view.
            typedef T_EXPR const& type;
          };
        //! Specialization for calls.
        template< class T_THIS, class T_EXPR >
          struct result< T_THIS(T_EXPR) >
          {
            //! Returns a view.
            typedef T_EXPR const& type;
          };
          
        //! implementation of the transform.
        template<typename T_EXPR >
          typename result<protect(const T_EXPR&)> :: type 
            operator ()( const T_EXPR& _expr ) const
              { return _expr; }
      };
    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif
