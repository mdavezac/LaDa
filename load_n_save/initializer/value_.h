//
//  Version: $Id: value_.h 1189 2009-06-17 02:19:52Z davezac $
//

#ifndef _LADA_LOADNSAVE_INITIALIZER_VALUE_H_
#define _LADA_LOADNSAVE_INITIALIZER_VALUE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/utility/result_of.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/fusion/include/fold.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/static_assert.hpp>

#include "../tree/tree.h"

#include "../grammar/grammar.h"
#include "../transforms/getname.h"
#include "../transforms/gettag.h"
#include "../transforms/getaction.h"
#include "../transforms/getvalues.h"
#include "values.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace initializer
    {
      class Value_
      {
          //! Type of the xml-like tree.
          typedef tree::Section t_Tree;
     
        public:
          //! Constructor.
          Value_( t_Tree::t_String const& _secname ) : secname_(_secname) {}
          //! Copy Constructor.
          Value_( Value_ const &_c ) : secname_(_c.secname_) {}
         
          //! Result type of the functor.
          typedef bool result_type;
     
          //! Functor.
          template< class  T_EXPR, class T_ACTION >
           bool operator()( T_EXPR const& _expr, T_ACTION& _action,
                            t_Tree::t_String const &_opname,
                            t_Tree::t_String const &_opvalue ) const;
          
        protected:
          //! Performs check and parsin.
          template< class T_ACTION > struct Check;
          //! Name of the section.
          t_Tree::t_String const& secname_;
      };

      template< class T_ACTION >
        struct Value_::Check : protected StraightValue_
        {
          typedef void result_type;
          Check   ( t_Tree::t_String const& _op, T_ACTION &_ac ) 
                : option_(_op), action_(_ac), found( false ), many(false) {}
          Check   ( Check const &_c )
                : option_(_c.option), action_(_c.action_), found(_c.found), many(_c.many) {}
          template< class  T_VALUE >
            void operator()( T_VALUE const &_value ) const
            {
              if( not StraightValue_()( option_, _value, action_ ) ) return;
              if( found ) many = true;
              found = true;
            }
  
          t_Tree::t_String const &option_;
          T_ACTION &action_;
          mutable bool found;
          mutable bool many;
        };
     
      template< class  T_EXPR, class T_ACTION >
        bool Value_::operator()( T_EXPR const& _expr, T_ACTION &_action,
                                 t_Tree::t_String const &_opname,
                                 t_Tree::t_String const &_opvalue ) const
        {
          namespace proto = boost::proto;

          Check< T_ACTION > check( _opvalue, _action );
          fusion::for_each( _expr, check );
          if( not check.found ) return false;
          if( check.many ) 
          {
            std::cerr << "More than one match for value " << _opvalue 
                      << " in option " << _opname << " of section " << secname_ << "\n"
                      << "Could be a bug in code.\n";
            return false;
          }
        
          return true;
        }
      
    } // namespace initializer.
  } // namespace load_n_save

} // namespace LaDa


#endif
