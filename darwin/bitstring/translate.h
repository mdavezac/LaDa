//
//  Version: $Id$
//
#ifndef _LADA_GA_BITSTRING_TRANSLATE_H_
#define _LADA_GA_BITSTRING_TRANSLATE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <ostream>

#include "print/manip.h"

namespace LaDa
{
  namespace GA
  {
    namespace BitString
    {
      //! Policy to translate an object containing a vector of 1:0 values to a structure.
      template< class T_OBJECT >
        struct Translate
        {
          //! Type of the Object.
          typedef T_OBJECT t_Object;
          //! From Crystal::Structure to objects.
          static void translate( const t_Object&, std :: string& );
          //! From Crystal::Structure to objects.
          static void translate( const std::string&, t_Object& );
        };

      template< class T_OBJECT >
        void Translate<T_OBJECT> :: translate( const std::string& _string, 
                                               t_Object& _object ) 
        {
          typedef typename t_Object :: t_Container :: iterator t_ivar;
          std::istringstream sstr( Print::StripEdges(_string) );
          _object.Container().clear();
          do
          {
            typename t_Object :: t_Container :: value_type type;
            sstr >> type;
            _object.Container().push_back(type);
          }
          while( not sstr.eof() );
        }

      template< class T_OBJECT >
        void Translate<T_OBJECT> :: translate( const t_Object& _object,
                                               std::string &_string ) 
        {
          typedef typename t_Object :: t_Container :: const_iterator t_ivar;
          std::ostringstream sstr( _string );
          t_ivar i_var = _object.Container().begin();
          t_ivar i_var_end = _object.Container().end();
          for(; i_var != i_var_end; ++i_var ) sstr << (*i_var) << " ";
          _string = Print::StripEdges(sstr.str()); 
        }

    } // namespace GroundStates.
  } // namespace GA
} // namespace LaDa 
#endif
