//
//  Version: $Id$
//
#ifndef _OPT_FRPRMN_H_
#define _OPT_FRPRMN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <nag.h>
#include <nage04.h>

#include "opt/opt_minimize_base.h"

//! \brief Type of the evaluation function an Minimizer::Nag will call during
//!        minimization.
typedef void (*t_NagFunction) ( int double*, double*, double*, Nag_Comm*);

//! \brief returns a pointer to the correct extern "C" evaluation function
//!        for Minimizer::Nag.
//! \details This routine should be specifically specialized for each
//!          functional which wants to use Minimizer::Nag
//! \sa Minimizer::typical_nagfun(), Vff::vff_nagfun(), Vff::layeredvff_nagfun().
template< class T_FUNCTIONAL > t_NagFunction choose_nag_function();

namespace LaDa
{
  namespace Minimizer {

    //! Should implement nag conjugate gradient. Don't have nag c hearers to try it!
    template<typename T_FUNCTIONAL>
    class Nag : public Base< T_FUNCTIONAL >
    {
      public:
        //! the type of the functional to minimize
        typedef T_FUNCTIONAL t_Functional;

        //! static member which c functions will call
        static t_Functional* static_pointer;

      protected:
        using Base<T_FUNCTIONAL>::current_func;
        int x_length; //!< length of x and grad ==  t_Functional :: variables->size()
        double tolerance; //!< Convergence for overall algorithm
        double line_tolerance; //!< for linmin
        int itermax; //!< maximum number of iterations
        //! Pointer to the evaluation function
        t_NagFunction func_call;
        //! How much to print
        types::t_unsigned doprint;
                                    
        
      public:
        //! Constructor 
        Nag () : x_length(current_func->size()),
                 tolerance(types::tolerance),
                 line_tolerance(types::tolerance),
                 itermax(500),
                 func_call(NULL), 
                 doprint(0) {}
        //! Constructor and Initializer
        Nag    (t_Functional& _func, t_NagFunction _f) 
             : Base<t_Functional>( _func ),
               x_length(0),
               tolerance(types::tolerance),
               line_tolerance(types::tolerance),
               zeps(types::tolerance),
               itermax(500),
               func_call(_f),
               doprint(0) {}
        //! Copy Constructor
        Nag   (const Nag& _c)
             : Base<t_Functional>( _c ),
               x_length(_c.x_length),
               tolerance(_c.tolerance),
               line_tolerance(_c.line_tolerance),
               itermax(_c.itermax),
               func_call(_c.func_call),
               doprint(_c.doprint) {}

        //! Destructor
        virtual ~Nag(){};

        //! Calls minimizer
        virtual bool operator()();

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
                             types::t_real _linetolerance );
        {
          itermax = _itermax;
          tolerance = _tolerance;
          line_tolerance = _linetolerance;
        }
    };

    //! \brief Transform an object into a "C" function callable by NAG
    //! \sa choose_nag_function(), Vff::vff_nagfun(), Vff::layeredvff_nagfun().
    template< class T_FUNCTIONAL >
    void typical_nagfun(int _n, double* _x, double* _r, double* _g, Nag_Comm*)
    {
      typedef typename Nag<T_FUNCTIONAL> t_Minimizer;

      if ( not (    t_Minimizer :: static_pointer 
                 or _r or _g or _x                ) ) return;

      const double *i_x_copy = _x;
      const double *i_x_end = _x + t_Minimizer::static_pointer->size();
      std::copy(i_x_copy, i_x_end, t_Minimizer::static_pointer->begin());
      *_r = t_Minimizer::static_pointer->evaluate_with_gradient(_g);
    }

    template<typename T_FUNCTIONAL> 
      T_FUNCTIONAL* Nag<T_FUNCTIONAL> :: static_pointer = NULL;
   

    template<typename T_FUNCTIONAL> 
    inline bool Nag<T_FUNCTIONAL> :: operator()()
    {
      // initializes object related stuff
      if( not current_func->init() ) return false;
      x_length = current_func->size();
      
      double *x, *g;
      try
      {
        x = new double[x_length];
        g = new double[x_length];
      }
      catch ( std::exception &_e )
      {
        std::ostringstream sstr;
        sstr << __FILE__ << ", line: " << __LINE__ << "\n"
             << "Could not allocate memory prior to calling Nag\n"
             << _e.what() << "\n";
        if( x ) delete[] x;
        if( g ) delete[] g;
        throw std::runtime_error( _e );
      }
      double* x_copy = x; 
      double result;

      // Sets up nag options
      static NagError fail;
      Nag_E04_Opt nagopts;
      e04xxc( &nagopt );

      nagopt.verify_grad    =  Nag_NoCheck;
      nagopt.max_iter       = (int) itermax;
      nagopt.optim_tol      = (double) tolerance;
      nagopt.linesearch_tol = (double) line_tolerance;
      switch( doprint )
      {
        case 0: nagopt.print_level = Nag_NoPrint; break;
        case 1: nagopt.print_level = Nag_Soln; break;
        case 2: nagopt.print_level = Nag_Iter; break;
        case 3: nagopt.print_level = Nag_Soln_Iter; break;
        default: nagopt.print_level = Nag_NoPrint; break;
      }

      std::copy ( current_func->begin(), current_func->end(), x_copy );
      static_pointer = current_func;
      e04dgc( x_length, nagfunc, result, g, &nagopts, NAGCOM_NULL, &fails);

      delete[] x;
      return fail.code == NE_NOERROR;
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
      _element.Attribute( "print", &doprint );

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
        if ( str.compare("nag" ) == 0 )
          break;
        parent = parent->NextSiblingElement("Minimizer");
      }
      if ( not parent )
      {
        std::cerr << "Could not find an <Minimizer type=\"frprmn\"> tag in input file" 
                  << std::endl;
        return false;
      } 
      return parent;
    }

    template<typename T_FUNCTIONAL> 
    inline bool Frpr<T_FUNCTIONAL> :: Load( const TiXmlElement &_node )
    {
      const TiXmlElement* parent = find_node( _node );
      if( parent ) return Load_(*parent);
      std::cerr << "Could not find an <Minimizer type=\"some gsl\"> tag in input file" 
                << std::endl;
      return false;
    }
  }
} // namespace LaDa
#endif
