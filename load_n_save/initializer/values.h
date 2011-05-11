//
//  Version: $Id: values.h 1226 2009-07-13 06:28:01Z davezac $
//

#ifndef _LADA_LOADNSAVE_INITIALIZER_VALUES_H_
#define _LADA_LOADNSAVE_INITIALIZER_VALUES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/xpressive/regex_compiler.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/joint_view.hpp>
#include <boost/fusion/include/as_vector.hpp>

#include "../tree/tree.h"

#include "../grammar/grammar.h"
#include "parse_value.h"
#include "convert_value_to_regex.h"
//#include "convert_value_to_action.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace initializer
    {
      namespace fusion = boost::fusion;

      class StraightValue_
      {
          //! Type of the xml-like tree.
          typedef tree::Section t_Tree;
     
        public:
          //! Constructor.
          StraightValue_() {}
         
          //! Result type of the functor.
          typedef bool result_type;
     
          //! Functor when neither value nor actoin exist.
          bool operator()( t_Tree::t_String const &_treeval,
                           boost::mpl::bool_<false> const&,
                           boost::mpl::bool_<false> const& ) const  { return _treeval.empty(); }
          //! Functor when value exists and is an std::string
          template< class T_ACTION >
            bool operator()( t_Tree::t_String const &_treeval,
                             std::string const &_exprval, 
                             T_ACTION &_action ) const
            {
              namespace bx = boost::xpressive;
              bx::sregex const re = bx::sregex::compile( std::string( _exprval ) );
              if(not bx::regex_match( _treeval, re ) ) return false;
              return parse_value( _treeval, _action );
            }
          //! Functor when value does not exist.
          template< class T_ACTION >
            bool operator()( t_Tree::t_String const &_treeval,
                             boost::mpl::bool_<false> const &, 
                             T_ACTION &_action ) const
              { return parse_value( _treeval, _action ); }
          //! Functor when value exists.
          template< class  T_VALUE, class T_ACTION >
            bool operator()( t_Tree::t_String const &_treeval,
                             T_VALUE const &_exprval, 
                             T_ACTION &_action ) const
            {
              std::string const val( convert_value_to_regex(_exprval) );
              boost::mpl::bool_<false> f;
              if( not (*this)( _treeval, val, f ) ) return false;
              return parse_value( _treeval, _action );
            }

          //! Functor for something == other value.
          template< class  T_VAL1, class T_VAL2, class T_ACTION >
            bool operator()( t_Tree::t_String const &_treeval,
                             boost::fusion::joint_view<T_VAL1, T_VAL2> const &_exprval, 
                             T_ACTION &_action ) const
            {
              namespace bf = boost::fusion;
              // First tries item 0.
              typedef typename bf::joint_view<T_VAL1, T_VAL2> t_joint;
              typedef typename bf::result_of::as_vector<t_joint> :: type t_vector;
              t_vector const exprval( bf::as_vector(_exprval) );
              std::string const val1( convert_value_to_regex( bf::at_c<0>(exprval) ) );
              std::string const val2( convert_value_to_regex( bf::at_c<1>(exprval) ) );
              boost::mpl::bool_<false> f;
              bool const parsed_val1( (*this)( _treeval, val1, f ) );
              bool const parsed_val2( (*this)( _treeval, val2, f ) );
              // It is implied that the second value can be converted by an action.
              if( not (parsed_val1 or parsed_val2) ) return false;
              return parse_value( val2, _action );
            }
         
          
          //! Do nothing
          template< class  T_VALUE, class T_ACTION >
            bool check( t_Tree::t_String const &_treeval,
                        T_VALUE const &_exprval, 
                        T_ACTION &_action ) const
            {
              T_ACTION a( _action );
              return (*this)( _treeval, _exprval, a );
            }
      };

    } // namespace initializer.
  } // namespace load_n_save

} // namespace LaDa


#endif
