//
//  Version: $Id$
//
#ifndef _GSL_SIMPLEX_H_
#define _GSL_SIMPLEX_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iomanip>
#include <gsl/gsl_multimin.h>
#include <functional>
#include <algorithm>
#include <tinyxml/tinyxml.h>
#include <boost/lambda/bind.hpp>

#include "algorithms.h"
#include "gsl.h"
#include "gsl_mins.h"

#ifdef _MPI
#  include <mpi/mpi_object.h>
#endif

namespace Minimizer {

  //! \cond
  namespace details
  {
    template< class T_FUNCTION >
      double gsl_function( const gsl_vector* _x, void* _data )
      {
        T_FUNCTION *_this = (T_FUNCTION*) _data;
        types::t_real result = (*_this)( _x->data );
        return result;
      }
  }
  //! \endcond
  //! \brief Minimizer interfaces for the Gnu Scientific Library
  //! \details Interface to the following algorithms:
  //!         - Fletcher-Reeves conjugate gradient
  //!         - Polak-Ribiere conjugate gradient 
  //!         - vector Broyden-Fletcher-Goldfarb-Shanno algorithm
  //!         - vector Broyden-Fletcher-Goldfarb-Shanno algorithm.
  //!           Second implementation, recommended by GSL manual...
  //!         - Steepest descent
  //!         .
  //! \xmlinput see TagMinimizer
  class Simplex 
  {
    public:
      types::t_real tolerance; //!< Complete convergence
      types::t_real stepsize; //!< line step
      types::t_unsigned itermax; //!< maximum number of iterations
      bool verbose; //!< Wether to print out during minimization
                                  
      
    public:
      //! Constructor and Initializer
      Simplex   ( types::t_unsigned _itermax,
                  types::t_real _tol, 
                  types::t_real _stepsize ) 
              : tolerance(_tol), 
                stepsize(_stepsize), itermax(_itermax),
                verbose(false) {}
            
      //! Constructor
      Simplex () : tolerance(types::tolerance),
                   stepsize(0.1), 
                   itermax(500), verbose(false) {}
      //! Destructor
      virtual ~Simplex(){};

      //! Non-XML way to set-up the minimizers.
      void set_parameters( types::t_unsigned _itermax,
                           types::t_real _tol, 
                           types::t_real _stepsize );
      //! Minimization functor
      template< class T_FUNCTION >
        typename T_FUNCTION :: t_Return
          operator()( const T_FUNCTION &_func,
                      typename T_FUNCTION :: t_Arg &_arg ) const;
      //! \brief Finds the node - if it is there - which describes this minimizer
      //! \details Looks for a \<Minimizer\> tag first as \a _node, then as a
      //!          child of \a _node. Different minimizer, defined by the
      //!          attribute types are allowed:
      const TiXmlElement* find_node( const TiXmlElement &_node );
      //! \brief Loads Minimizer directly from \a _node.
      //! \details If \a _node is not the correct node, the results are undefined.
      bool Load_( const TiXmlElement &_node );
      //! Loads the minimizer from XML
      bool Load( const TiXmlElement &_node );
      //! Serializes a structure.
      template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version);
  };

  template<typename T_FUNCTION> 
    typename T_FUNCTION :: t_Return 
      Simplex :: operator()( const T_FUNCTION &_func,
                             typename T_FUNCTION :: t_Arg &_arg ) const
    {
      namespace bl = boost::lambda;
      gsl_multimin_fminimizer *solver;
      gsl_vector *ss;

      __DEBUGTRYBEGIN
 
        if ( verbose ) std::cout << "Starting GSL minimization\n";
 
        gsl_multimin_function gsl_func;
        gsl_func.f = &details::gsl_function<const T_FUNCTION>;
        gsl_func.n = _arg.size();
        gsl_func.params = (void*) &_func;
        
        solver = gsl_multimin_fminimizer_alloc( gsl_multimin_fminimizer_nmsimplex,
                                                gsl_func.n);
        if (not solver) return false;
        
        ss = gsl_vector_alloc ( gsl_func.n );
        gsl_vector_set_all (ss, stepsize);
        ::Gsl::Vector x( _arg );
        gsl_multimin_fminimizer_set( solver, &gsl_func, (gsl_vector*) x, ss );
        
        int status;
        double newe(0), olde(0);
        types::t_unsigned iter( 0 );

 
        do
        {
          iter++;
          olde=newe;
          status = gsl_multimin_fminimizer_iterate (solver);
 
          if( status ) break;
 
          types::t_real size = gsl_multimin_fminimizer_size (solver);
          status = gsl_multimin_test_size (size, 1e-2);
          if( status == GSL_SUCCESS )
          {
            if( verbose )
              std::cout << "break on simplex small" << std::endl; 
            break; 
          }


          newe = gsl_multimin_fminimizer_minimum(solver);
          if( verbose )
            std::cout << "  Simplex Iteration " << iter
                      << ": minimum = " << newe
                      << "   size =  " << size << "\n";
        }
        while ( status == GSL_CONTINUE and ( itermax == 0 or iter < itermax ) );
        if( verbose and ( status != GSL_SUCCESS and iter != itermax ) ) 
          std::cout << "Error while minimizing with gsl: "
                    << gsl_strerror( status ) << ".\n";
 
        newe = gsl_multimin_fminimizer_minimum( solver );
        if ( verbose )
          std::cout << "Final Iteration: " << newe << std::endl;
    
        gsl_vector *minx = gsl_multimin_fminimizer_x( solver );
        types::t_unsigned i(0);
        opt::concurrent_loop
        (
          _arg.begin(), _arg.end(), i,
          bl::_1 = bl::bind( &gsl_vector_get, minx, bl::_2 )
        );
 
        gsl_multimin_fminimizer_free (solver);
        gsl_vector_free(ss);
 
        return newe;
 
      __DEBUGTRYEND
      (
        gsl_multimin_fminimizer_free (solver);
        gsl_vector_free(ss);,
        "Error encountered while minimizing with the GSL library\n"
      )
 
    }  // dummy minimizer


  template<class ARCHIVE>
    void Simplex :: serialize(ARCHIVE & _ar, const unsigned int _version)
    {
       _ar & itermax;
       _ar & tolerance;
       _ar & stepsize;
    }
  
}

#endif
