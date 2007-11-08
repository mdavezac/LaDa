//
//  Version: $Id$
//
#ifndef _OPT_MINIMIZE_BASE_H_
#define _OPT_MINIMIZE_BASE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

//! \brief Should contain all minimizers
//! \details Minimizers are expected to be accessible through a specific
//!          interface, and to have access to the functional to minimize
//!          through a specific interface as well. More explicitely, minimizers
//!          should able to load their parameters from XML. They should take
//!          the type of functional to minimize as a template. They should not
//!          expect any behavior from a functional which is not explicitely
//!          defined by funtion::Base. 
namespace minimizer {

  //! Declares the required behaviors for a minimizer.
  template<class T_FUNCTIONAL> 
  class Base
  {
    public:
      //! The type of functional to minimize.
      typedef T_FUNCTIONAL t_Functional;
    protected:
      //! A pointer to the functional to minimize
      t_Functional* current_func;

    public:
      //! Constructor
      Base() { current_func = NULL; }
      //! Load the parameters from XML
      Base(const TiXmlElement &_node ) { Load(_node); }
      //! Constructor and Initializer
      Base(t_Functional &_func) : current_func(&_func) {}
      //! Copy Constructor
      Base( Base<t_Functional> &_func) : current_func(&_func) {}
      //! Destructor
      virtual ~Base(){}

      // Should assign Base::current_func and perform a minimization.
      inline bool operator()( t_Functional &_func )
      { 
        current_func = &_func;
        return operator()();
      }
      //! Should perform the minimization
      bool operator()() { return true; }  
      //! Should load parameters from XML.
      virtual bool Load( const TiXmlElement &_node) { return true; }
      
      //! Sets the functional to minimize.
      void set_object( t_Functional &_func ) { current_func = &_func; }
  };


}
#endif 
