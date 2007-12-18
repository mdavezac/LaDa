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

#ifdef _DOXY_HOODWINKER_
  //! \ingroup Fortran
  //! \brief Original minimizer from nanopse
  //! \details Try and look for nanopse documentation if you can find it...
  extern "C" double frprmn ( double *, int *, double *,
                             double *, double *, int *,
                             double *, double *, int *,
                             double * );
#endif

  /** \brief Interface to the fortran minimizer 
   *  \details Not easy to call generic c++ from fortran.
   * You must define a function similar to the one below which 
   * is explicitely none generic...
   * \code
      extern "C" {
        double call_it( double *_x,
                        double *_y)
        {
          static Vff::Functional* this_func;
          if ( not _y  )
          {
            this_func = (Vff::Functional*) _x;
            return 0;
          }
     
          if ( not this_func )
            return 0;
     
          std::copy(_x, _x + this_func->size(), this_func->end())
          return this_func->evaluate_with_gradient(_y);
        }
      }
   *  \endcode
   *  for the example above, 
   *  first initialize  
   *  \code
      FRPR_minimization< Vff::Functional > minimizer( (double*)vff_object, &call_it )
   *  \endcode
   *  then 
   *  \code
      call_it( &vff_object, NULL )
   *  \endcode
   *  finally
   *  \code
      minimizer.minimize() 
   *  \endcode
   *  remember to repeat FROM step 2 when changing from one vff_object to
   *  another  */
  template<typename T_FUNCTIONAL>
  class Frpr : public Base< T_FUNCTIONAL >
  {
      //\cond
      typedef double (*t_func_call) (const double* const, double* const);
      //\endcond

    public:
      //! Type of the functional to minimize;
      typedef T_FUNCTIONAL t_Functional;
    protected:
      using Base<T_FUNCTIONAL>::current_func;
      int x_length; //!< length of x and grad ==  t_Functional :: variables->size()
      double tolerance; //!< Convergence for overall algorithm
      double line_tolerance; //!< for linmin
      double zeps; //!< for linmin
      int iter;    //!< number of performed iterations
      int itermax; //!< maximum number of iterations
      double fret; //!< minimum value of the function
      double rtol; //!< achieved tolerance
      //! Pointer to the evaluation function
      t_func_call func_call;
                                  
      
    public:
      //! Constructor 
      Frpr () : x_length(current_func->size()),
                tolerance(types::tolerance),
                line_tolerance(types::tolerance),
                zeps(types::tolerance),
                itermax(500),
                func_call(NULL) {}
      //! Constructor and Initializer
      explicit
        inline Frpr   (t_Functional& _func, t_func_call _f) 
                    : Base<t_Functional>( _func ),
                      x_length(0),
                      tolerance(types::tolerance),
                      line_tolerance(types::tolerance),
                      zeps(types::tolerance),
                      itermax(500),
                      func_call(_f) {}
      //! Destructor
      virtual ~Frpr(){};

      //! Calls minimizer
      virtual bool operator()();

      //! Load minimizer parameters from XML
      bool Load( const TiXmlElement &_element );
  };

  //! \cond
  extern "C" double FC_FUNC(frprmn, FRPRMN)
                           ( double *, int *, double *,
                             double *, double *, int *,
                             double *, double *, int *,
                             double * );
  //! \endcond

  template<typename T_FUNCTIONAL> 
  inline bool Frpr<T_FUNCTIONAL> :: operator()()
  {
    // initializes object related stuff
    if( not current_func->init() ) return false;
    x_length = current_func->size();

    double* x = new double[x_length];
    double* x_copy = x; 
    double result;

    std::copy ( current_func->begin(), current_func->end(), x_copy );
    FC_FUNC(frprmn, frprmn) ( x, &x_length, &tolerance,
                              &line_tolerance, &zeps,
                              &iter, &result, (double*)func_call, 
                              &itermax, &rtol );

    delete[] x;
    return rtol <= tolerance;
  }


  template<typename T_FUNCTIONAL> 
  inline bool Frpr<T_FUNCTIONAL> :: Load( const TiXmlElement &_element )
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
      std::cerr << "Could not find an <Minimizer type=\"frprmn\"> tag in input file" 
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

}
#endif
