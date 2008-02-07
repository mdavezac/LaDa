//
//  Version: $Id$
//
#ifndef _OPT_FRPRMN_H_
#define _OPT_FRPRMN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _DEBUG
#include <stdexcept>
#endif

#include <functional>

#include "minimize_base.h"
#include "opt_minimize_base.h"
#include "debug.h"

//! \brief Type of the evaluation function an Minimizer::Frpr will call during
//!        minimization.
typedef double (*t_FrprFunction) ( double*, double* );

//! \brief returns a pointer to the correct extern "C" evaluation function
//!        for Minimizer::Nag.
//! \details This routine should be specifically specialized for each
//!          functional which wants to use Minimizer::Nag
//! \sa Minimizer::typical_nagfun(), Vff::vff_nagfun(), Vff::layeredvff_nagfun().
template< class T_FUNCTIONAL > t_FrprFunction choose_frpr_function();

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
   *  \see vff_frprfun()
   * \xmlinput see TagMinimizer */
  template<typename T_FUNCTIONAL>
  class Frpr : public Base< T_FUNCTIONAL >
  {
    public:
      //! Type of the functional to minimize;
      typedef T_FUNCTIONAL t_Functional;
      //! static member which c functions will call
      static t_Functional* static_pointer;

    protected:
      typedef Base<t_Functional> t_Base;
      using t_Base::current_func;
      int x_length; //!< length of x and grad ==  t_Functional :: variables->size()
      double tolerance; //!< Convergence for overall algorithm
      double line_tolerance; //!< for linmin
      double zeps; //!< for linmin
      int iter;    //!< number of performed iterations
      int itermax; //!< maximum number of iterations
      double fret; //!< minimum value of the function
      double rtol; //!< achieved tolerance
      //! Pointer to the evaluation function
      t_FrprFunction func_call;
                                  
      
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
        inline Frpr   (t_Functional& _func, t_FrprFunction _f) 
                    : Base<t_Functional>( _func ),
                      x_length(0),
                      tolerance(types::tolerance),
                      line_tolerance(types::tolerance),
                      zeps(types::tolerance),
                      itermax(500),
                      func_call(_f) {}
      explicit
        inline Frpr   (t_FrprFunction _f) 
                    : Base<t_Functional>(),
                      x_length(0),
                      tolerance(types::tolerance),
                      line_tolerance(types::tolerance),
                      zeps(types::tolerance),
                      itermax(500),
                      func_call(_f) {}
      //! Copy Constructor
      Frpr   (const Frpr& _c)
           : Base<t_Functional>( _c ),
             x_length(_c.x_length),
             tolerance(_c.tolerance),
             line_tolerance(_c.line_tolerance),
             zeps(_c.zeps),
             itermax(_c.itermax),
             func_call(_c.func_call) {}

      //! Destructor
      virtual ~Frpr(){};

      //! Calls minimizer
      virtual bool operator()();
      using t_Base::operator();

      //! Load minimizer parameters from XML
      bool Load( const TiXmlElement &_element );

      //! \brief Finds the node - if it is there - which describes this minimizer
      //! \details Looks for a \<Minimizer\> tag first as \a _node, then as a
      //!          child of \a _node. Different minimizer, defined by the
      //!          attribute types are allowed:
      const TiXmlElement* find_node( const TiXmlElement &_node );

      //! \brief Loads Minimizer directly from \a _node.
      //! \details If \a _node is not the correct node, the results are undefined.
      bool Load_( const TiXmlElement &_node );

      //! sets all parameters
      void set_parameters( types::t_unsigned _itermax,
                           types::t_real _tolerance,
                           types::t_real _linetolerance,
                           types::t_real _zeps )
      {
        itermax = _itermax;
        tolerance = _tolerance;
        line_tolerance = _linetolerance;
        zeps = _zeps;
      }
  };
 
  //! \brief Transform an object into a "C" function callable by Minimizer::Frpr
  //! \sa choose_frpr_function(), Vff::vff_frprfun(), Vff::layeredvff_frprfun().
  template< class T_FUNCTIONAL >
  double typical_frprfun(double* _x, double* _y)
  {
    typedef Frpr<T_FUNCTIONAL> t_Minimizer;

    if ( not ( t_Minimizer :: static_pointer or _x or _y ) ) return 0;

    const double *i_x_copy = _x;
    const double *i_x_end = _x + t_Minimizer::static_pointer->size();
    std::copy(i_x_copy, i_x_end, t_Minimizer::static_pointer->begin());
    double result = t_Minimizer::static_pointer->evaluate_with_gradient(_y);
    return result;
  }


  template<typename T_FUNCTIONAL> 
    T_FUNCTIONAL* Frpr<T_FUNCTIONAL> :: static_pointer = NULL;
 

  //! \cond
  extern "C" double FC_FUNC(frprmn, FRPRMN)
                           ( double *, int *, double *,
                             double *, double *, int *,
                             double *, void *, int *,
                             double * );
  //! \endcond

  template<typename T_FUNCTIONAL> 
  inline bool Frpr<T_FUNCTIONAL> :: operator()()
  {
    // initializes object related stuff
    __ASSERT(not current_func,
             "Pointer to the function to minimize has not been set.\n")

    if( not current_func->init() ) return false;
    x_length = current_func->size();

    double* x;
    __TRYDEBUGCODE( x = new double[x_length];,
                   "Memory Allocation error.\n")

#define _MXPARM_ 300
    __DOASSERT( x_length > _MXPARM_,
                  "Current optimizer cannot go beyond " << _MXPARM_ << " variables\n"
               << "Change file cited above, as well as variable mxparm in opt/df1dim.f90 "
               << "and opt/linmin.f90 and recompile if you need to optimize larger structures.\n";)
    double* x_copy = x; 
    double result;

    static_pointer = current_func;
    std::copy ( current_func->begin(), current_func->end(), x_copy );
    FC_FUNC(frprmn, frprmn) ( x, &x_length, &tolerance,
                              &line_tolerance, &zeps,
                              &iter, &result, (void*)func_call, 
                              &itermax, &rtol );

    delete[] x;
    return rtol <= tolerance;
  }


  template<typename T_FUNCTIONAL> 
  inline bool Frpr<T_FUNCTIONAL> :: Load_( const TiXmlElement &_element )
  {
    _element.Attribute( "itermax", &itermax );
    if ( not itermax ) itermax = 500;
    _element.Attribute( "tolerance", &tolerance );
    if ( not tolerance ) tolerance = types::tolerance;
    _element.Attribute( "linetolerance", &line_tolerance);
    if ( not line_tolerance ) line_tolerance = tolerance;
    _element.Attribute( "zeps", &zeps );
    if ( not zeps ) zeps = tolerance;
  
    return true;
  }


  template<typename T_FUNCTIONAL> 
  inline const TiXmlElement* Frpr<T_FUNCTIONAL> :: find_node( const TiXmlElement &_node )
  {
    const TiXmlElement *parent;
    std::string str;
  
    // This whole section tries to find a <Functional type="vff"> tag
    // in _element or its child
    str = _node.Value();
    if ( str.compare("Minimizer" ) != 0 )
      parent = _node.FirstChildElement("Minimizer");
    else
      parent = &_node;
  
    
    while (parent)
    {
      str = "";
      if ( parent->Attribute( "type" )  )
        str = parent->Attribute("type");
      if ( str.compare("frprmn" ) == 0 )
        break;
      parent = parent->NextSiblingElement("Minimizer");
    }
    __DOASSERT( not parent, 
                "Could not find an <Minimizer type=\"frprmn\"> tag in input file.\n" )
    return parent;
  }

  template<typename T_FUNCTIONAL> 
  inline bool Frpr<T_FUNCTIONAL> :: Load( const TiXmlElement &_node )
  {
    const TiXmlElement* parent = find_node( _node );
    if( parent ) return Load_(*parent);
    return false;
  }
}
#endif
