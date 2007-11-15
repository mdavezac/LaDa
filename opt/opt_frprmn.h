//
//  Version: $Id$
//
#ifndef _OPT_FRPRMN_H_
#define _OPT_FRPRMN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <functional>

#include "opt/opt_minimize_base.h"

namespace Minimizer {

  extern "C" { double frprmn_( double *, int *, double *,
                               double *, double *, int *,
                               double *, double *, int *,
                               double * ); }

// Not easy to call generic c++ from fortran.
// You must define a function similar to the one below which 
// is explicitely none generic...
//
//
//  extern "C" {
//    double call_it( double *_x,
//                    double *_y)
//    {
//      static Vff::Functional* this_func;
//      if ( not _y  )
//      {
//        this_func = (Vff::Functional*) _x;
//        return 0;
//      }
// 
//      if ( not this_func )
//        return 0;
// 
//      std::copy(_x, _x + this_func->size(), this_func->end())
//      return this_func->evaluate_with_gradient(_y);
//    }
//  }
//
//  for the example above, 
//  first initialize  
//  FRPR_minimization< Vff::Functional > minimizer( (double*)vff_object, &call_it )
//  then 
//  call_it( &vff_object, NULL )
//  finally
//  minimizer.minimize() 
//  remember to repeat FROM step 2 when changing from one vff_object to
//  another 

  template<typename Object> 
  class Frpr : public Base< Object >
  {
      typedef double (*t_func_call) (const double* const, double* const);
    protected:
      using Base<Object>::current_func;
      int x_length; // length of x and grad ==  Object :: variables->size()
      double tolerance;
      double line_tolerance, zeps; // for linmin
      int iter; // number of performed iterations
      int itermax; // maximum number of iterations
      double fret; // minimum value of the function
      double rtol; // achieved tolerance
      t_func_call func_call;
                                  
      
    public:
      Frpr () : x_length(current_func->size()),
                tolerance(types::tolerance),
                line_tolerance(types::tolerance),
                zeps(types::tolerance),
                itermax(500),
                func_call(NULL) {}
      explicit
        inline Frpr   (Object& _func, t_func_call _f) 
                    : Base<Object>( _func ),
                      x_length(current_func->size()),
                      tolerance(types::tolerance),
                      line_tolerance(types::tolerance),
                      zeps(types::tolerance),
                      itermax(500),
                      func_call(_f) {}
      virtual ~Frpr(){};

      virtual bool minimize()
      {
        // initializes object related stuff
        if( not current_func->init() )
          return false;
        x_length = current_func->size();

        double* x = new double[x_length];
        double* x_copy = x; 
        double result;

        std::copy ( current_func->begin(), current_func->end(), x_copy );

        frprmn_( x, &x_length, &tolerance, &line_tolerance, &zeps,
                 &iter, &result, (double*)func_call, &itermax, &rtol );

        delete[] x;
        return rtol <= tolerance;
      }  // dummy minimizer


      bool Load( const TiXmlElement &_element )
      {
        const TiXmlElement *parent;
        std::string str;
      
        // This whole section tries to find a <Functional type="vff"> tag
        // in _element or its child
        str = _element.Value();
        if ( str.compare("Minimizer" ) != 0 )
          parent = _element.FirstChildElement("Minimizer");
        else
          parent = &_element;
      
        
        while (parent)
        {
          str = "";
          if ( parent->Attribute( "type" )  )
            str = parent->Attribute("type");
          if ( str.compare("frprmn" ) == 0 )
            break;
          parent = parent->NextSiblingElement("Minimizer");
        }
        if ( not parent )
        {
          std::cerr << "Could not find an <Minimizer type=\"FRP\"> tag in input file" 
                    << std::endl;
          return false;
        } 
      
        parent->Attribute( "itermax", &itermax );
        if ( not itermax ) itermax = 500;
        parent->Attribute( "tolerance", &tolerance );
        if ( not tolerance ) tolerance = types::tolerance;
        parent->Attribute( "linetolerance", &line_tolerance);
        if ( not line_tolerance ) line_tolerance = tolerance;
        parent->Attribute( "zeps", &zeps );
        if ( not zeps ) zeps = tolerance;
      
        return true;
     }
  };

}
#endif
