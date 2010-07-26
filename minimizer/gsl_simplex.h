#ifndef _GSL_SIMPLEX_H_
#define _GSL_SIMPLEX_H_

#include "LaDaConfig.h"

#include <iomanip>
#include <functional>
#include <algorithm>

#include <boost/tuple/tuple.hpp>
#include <boost/lambda/bind.hpp>

#include <tinyxml/tinyxml.h>
#include <gsl/gsl_multimin.h>

#include <opt/algorithms.h>

#include "gsl.h"
#include "gsl_mins.h"

#ifdef LADA_MPI
#  include <mpi/mpi_object.h>
#endif

namespace LaDa
{
  namespace Minimizer 
  {


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
        friend class boost::serialization::access;
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
                        typename T_FUNCTION :: t_Arg &_arg ) const
            { return operator_< T_FUNCTION,
                                typename T_FUNCTION :: t_Arg,
                                typename T_FUNCTION :: t_Return
                              >( _func, _arg ); }
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
      private:
        //! Minimization functor
        template< class T_FUNCTION, class T_CONTAINER, class T_RETURN >
          T_RETURN operator_( const T_FUNCTION &_func, T_CONTAINER &_arg ) const;
        //! Serializes a structure.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version);
    };

    LADA_REGISTER_MINIMIZER_VARIANT_HEADER( Simplex, "GSL Simplex" )

    template<class T_FUNCTION, class T_CONTAINER, class T_RETURN> 
      T_RETURN  Simplex :: operator_( const T_FUNCTION &_func, T_CONTAINER &_arg ) const
      {
        namespace bl = boost::lambda;
        gsl_multimin_fminimizer *solver;
        gsl_vector *ss;
        typedef boost::tuples::tuple< const T_FUNCTION&, T_CONTAINER& > t_Pair;
        t_Pair data_pair( _func, _arg );

        LADA_DEBUG_TRY_BEGIN
   
          if ( verbose ) std::cout << "Starting GSL minimization\n";
   
          gsl_multimin_function gsl_func;
          gsl_func.f = &details::gsl_f<t_Pair>;
          gsl_func.n = _arg.size();
          gsl_func.params = (void*) &data_pair;
          
          solver = gsl_multimin_fminimizer_alloc( gsl_multimin_fminimizer_nmsimplex,
                                                  gsl_func.n);
          if (not solver) return false;
          
          ss = gsl_vector_alloc ( gsl_func.n );
          gsl_vector_set_all (ss, stepsize);
          LaDa::Gsl::Vector x( _arg );
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
   
        LADA_DEBUG_TRY_END
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
} // namespace LaDa

#endif
