//
//  Version: $Id$
//
#ifndef _OPT_MINIMIZE_BASE_H_
#define _OPT_MINIMIZE_BASE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

namespace minimizer {

  template<class T_FUNCTIONAL> 
  class Base
  {
    public:
      typedef T_FUNCTIONAL t_Functional;
    // in order to interact with c codes, this object is static
    protected:
      t_Functional* current_func;

    public:
      inline Base() 
        { current_func = NULL; }
      inline Base(const TiXmlElement &_node )
        { Load(_node); }
      inline Base(t_Functional &_func) 
        { current_func = &_func; }
      virtual ~Base(){}

      // should contain at least two functions:
      inline bool operator()( t_Functional &_func )
      { 
        current_func = &_func;
        return operator()();
      }
      virtual bool operator()()
        { return true; }  // dummy minimizer
      virtual bool Load( const TiXmlElement &_node) = 0;
      
      void set_object( t_Functional &_func )
        { current_func = &_func; }
  };


}
#endif 
