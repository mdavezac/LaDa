//
//  Version: $Id$
//
#ifndef _DARWIN_GROUNDSTATES_TRANSLATE_H_
#define _DARWIN_GROUNDSTATES_TRANSLATE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lexical_cast.hpp>
#include <string>
#include <ostream>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <crystal/structure.h>
#include <opt/types.h>
#include <mpi/mpi_object.h>
#include <vff/layered.h>

#include "evaluator_base.h"
#include "../bitstring.h"
#include "../ce.h"

#include "translate.h"
#include "assign.h"

namespace LaDa
{
  namespace GA
  {
    namespace GroundStates
    {
      //! Policy to translate an object containing a vector of 1:0 values to a structure.
      template< class T_OBJECT >
        struct Translate
        {
          //! Type of the Object.
          typedef T_OBJECT t_Object;
          //! From objects to Crystal::Structure. 
          static void translate( const t_Object&, Crystal :: Structure& );
          //! From Crystal::Structure to objects.
          static void translate( const Crystal :: Structure&, t_Object& );
          //! From Crystal::Structure to objects.
          static void translate( const t_Object&, std :: string& );
          //! From Crystal::Structure to objects.
          static void translate( const std::string&, t_Object& );
        };

      template< class T_OBJECT >
        void Translate<T_OBJECT> :: translate( const Crystal::Structure& _structure,
                                               t_Object& _object )
        {
          typedef typename t_Object :: t_Container :: iterator t_ivar;
          container.clear();

          typedef typename Crystal::Structure::t_Atoms::const_iterator t_iatom;
          typedef typename t_Object :: t_Container :: iterator t_ivar;
          t_iatom i_atom = _structure.atoms.begin();
          t_iatom i_atom_end = _structure.atoms.end();
          for( ; i_atom != i_atom_end; ++i_atom, ++i_var )
            *i_var = i_atom->type > 0 ? 1: 0;
        }

      template< class T_OBJECT >
        void Translate<T_OBJECT> :: translate( const t_Object& _object, 
                                               Crystal::Structure& _structure )
        {
          typedef typename t_Object :: t_Container :: iterator t_ivar;
          typedef typename Crystal::Structure::t_Atoms::const_iterator t_iatom;
          typedef typename t_Object :: t_Container :: iterator t_ivar;
          t_iatom i_atom = _structure.atoms.begin();
          t_iatom i_atom_end = _structure.atoms.end();
          for( ; i_atom != i_atom_end; ++i_atom, ++i_var )
            i_atom->type = *i_var > 0 ? 1.e0: -1e0;
        }

      template< class T_OBJECT >
        void Translate<T_OBJECT> :: translate( const std::string& _string, 
                                               t_Object& _object ) 
        {
          typedef typename t_Object :: t_Container :: iterator t_ivar;
          std::istringstream sstr( _string );
          t_ivar i_var = _object.Container().begin();
          __DODEBUGCODE( t_ivar i_var_end = _object.Container().end(); )
          do
          {
            __ASSERT( i_var == i_var_end, "Object smaller than string.\n" )
            sstr >> *i_var;
            ++i_var;
          }
          while( not sstr.eof() );
          __ASSERT( i_var != i_var_end, "String smaller than object.\n" )
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
          _string = sstr.str(); 
        }

    } // namespace GroundStates.
  } // namespace GA
} // namespace LaDa 
#endif
