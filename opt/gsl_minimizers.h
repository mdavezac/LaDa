//
//  Version: $Id$
//
#ifndef _GSL_MINIMIZERS_H_
#define _GSL_MINIMIZERS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gsl/gsl_multimin.h>
#include <functional>
#include <algorithm>

#include "minimize_base.h"

#ifdef _MPI
#  include <mpi/mpi_object.h>
#endif

namespace LaDa
{
  namespace Minimizer {

    //! Creates an iterator which covers the gsl variables
    class gsl_iterator
    {
      public:
        //! type of the iterator
        typedef double* iterator_category;
        //! type of the value which the iterator references
        typedef double value_type;
        //! type of the difference of two iterators
        typedef int difference_type;
        //! type of the iterator
        typedef double* pointer;
        //! type of the value which the iterator references
        typedef double reference;

      private:
        //! The vector over which to iterate
        gsl_vector *x;
        //! Position in the vector
        double* iterator;

      public:
        //! Constructor and Initializer
        gsl_iterator( gsl_vector *_x ) : x(_x), iterator(x->data) {};
        //! Constructor and Initializer
        gsl_iterator( const gsl_iterator &_x ) : x(_x.x), iterator(_x.iterator) {};
        //! Moves iterator position forward by 1. (pre-incrementation)
        gsl_iterator& operator++()
         { iterator += x->stride; return *this; }
        //! Moves iterator position forward by 1. (post-incrementation)
        gsl_iterator& operator++(int)
         { iterator += x->stride; return *this; }
        //! Returns an iterator moved \a n positions forward.
        gsl_iterator operator+(types::t_unsigned n)
         { gsl_iterator it(x); it += n; return it;}
        //! Moves iterator position backward by \a n. 
        gsl_iterator& operator+=(types::t_unsigned n)
         { iterator += n*x->stride; return *this; }
        //! Returns an iterator moved \a n positions backward.
        gsl_iterator operator-(types::t_unsigned n)
         { gsl_iterator it(x); it -= n; return it;}
        //! Moves iterator position backward by 1. (pre-decrementation)
        gsl_iterator& operator--()
         { iterator -= x->stride; return *this; }
        //! Moves iterator position backward by \a n.
        gsl_iterator& operator-=(types::t_unsigned n)
         { iterator -= n*x->stride; return *this; }
        //! Returns a reference to the value at the current position
        double& operator*() const
         { return *iterator; }
        //! Returns true \a _x is at the same position
        bool operator==(const gsl_iterator &_x) const
         { return iterator == _x.iterator; }
        //! Returns true \a _x is not at the same position
        bool operator!=(const gsl_iterator &_x) const
         { return iterator != _x.iterator; }
    };
    //! Creates a constant iterator which covers the gsl variables
    class gsl_const_iterator
    {
      public:
        //! type of the constant iterator
        typedef const double* iterator_category;
        //! type of the constant value which the iterator references
        typedef const double value_type;
        //! type of the difference of two iterators
        typedef int difference_type;
        //! type of the constant iterator
        typedef const double* pointer;
        //! type of the value which the iterator references
        typedef const double reference;

      private:
        //! The constant vector over which to iterate
        const gsl_vector *x;
        //! Position in the constant vector
        pointer iterator;
      public:
        //! Constructor and Initializer
        gsl_const_iterator( const gsl_vector *_x ) : x(_x), iterator(x->data) {};
        //! Constructor and Initializer
        gsl_const_iterator   ( const gsl_const_iterator &_x ) 
                           : x(_x.x), iterator(_x.iterator) {};
        //! Moves constant iterator position forward by 1. (pre-incrementation)
        gsl_const_iterator& operator++()
         { iterator += x->stride; return *this; }
        //! Moves iterator position forward by 1. (post-incrementation)
        gsl_const_iterator& operator++(int)
         { iterator += x->stride; return *this; }
        //! Moves constant iterator position forward by \a n.
        gsl_const_iterator& operator+=(types::t_unsigned n)
         { iterator += n*x->stride; return *this; }
        //! Returns a constant iterator moved \a n positions forward.
        gsl_const_iterator operator+(types::t_unsigned n)
         { gsl_const_iterator it(x); it += n; return it;}
        //! Returns a constant iterator moved \a n positions backward.
        gsl_const_iterator operator-(types::t_unsigned n)
         { gsl_const_iterator it(x); it -= n; return it;}
        //! Moves constant iterator position backward by 1. (pre-decrementation)
        gsl_const_iterator& operator--()
         { iterator -= x->stride; return *this; }
        //! Moves constant iterator position backward by 1. (post-decrementation)
        gsl_const_iterator& operator--(int)
         { iterator -= x->stride; return *this; }
        //! Moves constant iterator position backward by \a n.
        gsl_const_iterator& operator-=(types::t_unsigned n)
         { iterator -= n*x->stride; return *this; }
        //! Returns a reference to the constant value at the current position
        reference& operator*() const
         { return *iterator; }
        //! Returns true \a _x is at the same position
        bool operator==(const gsl_const_iterator &_x)
         { return iterator == _x.iterator; }
        //! Returns true \a _x is not at the same position
        bool operator!=(const gsl_const_iterator &_x)
         { return iterator != _x.iterator; }
    };

    //! Interface for evaluating the functioanl
    template<typename T_FUNTIONAL>
    double gsl_func_f( const gsl_vector *_x, void *_params);
    //! Interface for evaluating the derivatives of functional
    template<typename T_FUNTIONAL>
    void gsl_func_df( const gsl_vector *_x, void *_params, gsl_vector *_grad);
    //! Interface for evaluating the functional and its derivatives.
    template<typename T_FUNTIONAL>
    void gsl_func_fdf( const gsl_vector *_x, void *_params, double *_r, gsl_vector *_grad);


    //! \brief Minimizer interfaces for the Gnu Scientific Library
    //! \details Interface to the following algorithms:
    //!         - Fletcher-Reeves conjugate gradient
    //!         - Polak-Ribiere conjugate gradient 
    //!         - vector Broyden-Fletcher-Goldfarb-Shanno algorithm
    //!         - vector Broyden-Fletcher-Goldfarb-Shanno algorithm.
    //!           Second implementation, recommended by GSL manual...
    //!         - Steepest descent
    //!         .
    //! \todo What about the simplex method...
    //! \xmlinput see TagMinimizer
    template<typename T_FUNCTIONAL> 
    class GnuSL : public Base< T_FUNCTIONAL >
    {
        friend class boost::serialization::access;
      public:
        typedef T_FUNCTIONAL t_Functional;


      protected:
        //!< Lists all known gsl multidimensional minimizers.
        enum t_gsl_minimizer_type
        { 
          GSL_NONE,   //!< No minizer... 
          GSL_FR,     //!< Fletcher-Reeves conjugate gradient algorithm.
          GSL_PR,     //!< Polak-Ribiere conjugate gradient algorithm.                  
          GSL_BFGS2,  //!< More efficient Broyden-Fletcher-Goldfarb-Shanno algorithm.
          GSL_BFGS,   //!< Broyden-Fletcher-Goldfarb-Shanno algorithm.
          GSL_SD      //!< Steepest Descent algorithm.
        };

      public:
        //! Fletcher-Reeves conjugate gradient algorithm.
        const static t_gsl_minimizer_type FletcherReeves = GSL_FR;
        //! Polak-Ribiere conjugate gradient algorithm.                  
        const static t_gsl_minimizer_type PolakRibiere = GSL_PR;
        //!  More efficient Broyden-Fletcher-Goldfarb-Shanno algorithm.
        const static t_gsl_minimizer_type BFGS2 = GSL_BFGS2;
        //!  Broyden-Fletcher-Goldfarb-Shanno algorithm.
        const static t_gsl_minimizer_type BFGS = GSL_BFGS;
        //!  Steepest Descent algorithm.
        const static t_gsl_minimizer_type SteepestDescent = GSL_SD;

      protected:
        typedef Base<t_Functional> t_Base;
        using t_Base::current_func;
        types::t_real tolerance; //!< Complete convergence
        types::t_real linetolerance; //!< Line convergences
        types::t_real linestep; //!< line step
        types::t_unsigned iter; //!< number of performed iterations
        types::t_unsigned itermax; //!< maximum number of iterations
        t_gsl_minimizer_type type; //!< maximum number of iterations
      public:
        bool do_print; //!< Wether to print out during minimization
                                    
        
      public:
        //! Constructor and Initializer
        GnuSL   ( t_gsl_minimizer_type _type, 
                  types::t_unsigned _itermax,
                  types::t_real _tol, 
                  types::t_real _linetol, 
                  types::t_real _linestep ) 
              : tolerance(_tol), linetolerance(_linetol),
                linestep(_linestep), iter(0), itermax(_itermax),
                type( _type ), do_print(false) {}
              
        //! Constructor
        GnuSL () : tolerance(types::tolerance),
                   linetolerance(0.01),
                   linestep(0.1), iter(0),
                   itermax(500), do_print(false) {}
        //! Constructor and Initializer
        explicit
          inline GnuSL   (t_Functional& _func)
                       : Base<t_Functional>( _func ),
                         tolerance(types::tolerance),
                         linetolerance(0.01),
                         linestep(0.1), iter(0),
                         itermax(500), do_print(false) {}
        //! Destructor
        virtual ~GnuSL(){};

        //! Non-XML way to set-up the minimizers.
        void set_parameters( t_gsl_minimizer_type _type, 
                             types::t_unsigned _itermax,
                             types::t_real _tol, 
                             types::t_real _linetol, 
                             types::t_real _linestep );
        //! Minimization functor
        virtual bool operator()();
        using t_Base::operator();
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
        //! Serializes a structure.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version);
    };

    template<typename T_FUNCTIONAL> 
    inline void GnuSL<T_FUNCTIONAL> :: set_parameters( t_gsl_minimizer_type _type, 
                                                       types::t_unsigned _itermax,
                                                       types::t_real _tol, 
                                                       types::t_real _linetol, 
                                                       types::t_real _linestep ) 
    {
      tolerance     = _tol;
      linetolerance = _linetol;
      linestep      = _linestep;
      iter          = 0;
      itermax       = _itermax;
      type          = _type;
    }
              

    template<typename T_FUNCTIONAL> 
    bool GnuSL<T_FUNCTIONAL> :: operator()()
    {
      const gsl_multimin_fdfminimizer_type *T;
      gsl_multimin_fdfminimizer *s;
      gsl_vector *x;
      int status;

      __DEBUGTRYBEGIN

        // initializes object related stuff
        if( not current_func )  return false;
        if( not current_func->init() )  return false;
        
        gsl_multimin_function_fdf gsl_func;
        
        gsl_func.f = &gsl_func_f<t_Functional>;
        gsl_func.df = &gsl_func_df<t_Functional>;
        gsl_func.fdf =  &gsl_func_fdf<t_Functional>;
        gsl_func.n = current_func->size();
        gsl_func.params = (void*) current_func;
        
        x = gsl_vector_alloc (gsl_func.n);
        gsl_iterator i_x(x);
        std::copy(current_func->begin(), current_func->end(), i_x );
        
        switch( type )
        {
          default:
          case GSL_BFGS2: T = gsl_multimin_fdfminimizer_vector_bfgs2; break;
          case GSL_FR: T = gsl_multimin_fdfminimizer_conjugate_fr; break;
          case GSL_PR: T = gsl_multimin_fdfminimizer_conjugate_pr; break;
          case GSL_BFGS: T = gsl_multimin_fdfminimizer_vector_bfgs; break;
          case GSL_SD: T = gsl_multimin_fdfminimizer_steepest_descent; break;
        }
        s = gsl_multimin_fdfminimizer_alloc (T, gsl_func.n);
        if (not s) return false;
        
        gsl_multimin_fdfminimizer_set (s, &gsl_func, x, linestep, linetolerance);
        
        double newe=current_func->energy(), olde=0;
        types::t_unsigned nomove=0;
        if ( do_print ) std::cout << "Starting GSL minimization\n";
        iter = 0;
      
        do
        {
          iter++;
          olde=newe;
          status = gsl_multimin_fdfminimizer_iterate (s);
        
          if (status)
          {
            printf ("error: %s\n", gsl_strerror (status));
            break;
          }
        
          newe = gsl_multimin_fdfminimizer_minimum(s);
          if ( tolerance > std::abs(newe - olde) )  ++nomove;
          else                                      nomove = 0;
          if ( do_print )
            std::cout << "Iteration " << iter 
                      << ": " << current_func->energy() << "\n";
        
        }
        while ( nomove < 3 and iter < itermax );

        gsl_func_f<t_Functional>( gsl_multimin_fdfminimizer_x(s), (void*)current_func );
        if ( do_print )
          std::cout << "Final Iteration: " << current_func->energy() << "\n\n";
    
        gsl_multimin_fdfminimizer_free (s);
        gsl_vector_free (x);

      __DEBUGTRYEND
      (
        gsl_multimin_fdfminimizer_free (s);
        gsl_vector_free (x);,
        "Error encountered while minimizing with the GSL library\n"
      )

      return status == GSL_SUCCESS;
    }  // dummy minimizer


    template<typename T_FUNCTIONAL> 
    const TiXmlElement* GnuSL<T_FUNCTIONAL> :: find_node( const TiXmlElement &_node )
    {
      const TiXmlElement *parent;
      std::string str;
    
      // This whole section tries to find a <Functional type="vff"> tag
      // in _node or its child
      str = _node.Value();
      if ( str.compare("Minimizer" ) != 0 )
        parent = _node.FirstChildElement("Minimizer");
      else parent = &_node;
    
      
      while (parent)
      {
        str = "";
        type = GSL_NONE;
        if ( parent->Attribute( "type" )  )
          str = parent->Attribute("type");
        if ( str.compare("gsl_bfgs2" ) == 0 )     type = GSL_BFGS2;
        else if ( str.compare("gsl_fr" ) == 0 )   type = GSL_FR;
        else if ( str.compare("gsl_pr" ) == 0 )   type = GSL_PR;
        else if ( str.compare("gsl_bfgs" ) == 0 ) type = GSL_BFGS;
        else if ( str.compare("gsl_sd" ) == 0 )   type = GSL_SD;
        if ( type != GSL_NONE ) break;
        parent = parent->NextSiblingElement("Minimizer");
      }
      if ( parent ) return parent;
      
      return NULL;
    }
    
    template<typename T_FUNCTIONAL> 
    bool GnuSL<T_FUNCTIONAL> :: Load_( const TiXmlElement &_node )
    {
      double d; int n;
      _node.Attribute( "itermax", &n );
      itermax = (n > 0) ? abs(n) : 10000;
      _node.Attribute( "tolerance", &d );
      tolerance = d ? types::t_real(d) : types::tolerance;
      _node.Attribute( "linetolerance", &d );
      linetolerance = d ? types::t_real(d) : 0.01;
      _node.Attribute( "linestep", &d );
      linestep = d ? types::t_real(d) : 0.1;
      if( _node.Attribute("doprint") ) do_print = true;
    
      return true;
    }
    
    template<typename T_FUNCTIONAL> 
    bool GnuSL<T_FUNCTIONAL> :: Load( const TiXmlElement &_node )
    {
      const TiXmlElement* parent = find_node( _node );
      if( parent ) return Load_(*parent);
      std::cerr << "Could not find an <Minimizer type=\"some gsl\"> tag in input file" 
                << std::endl;
      return false;
    }

    template<typename T_FUNCTIONAL> template<class ARCHIVE>
      void GnuSL<T_FUNCTIONAL> :: serialize(ARCHIVE & _ar, const unsigned int _version)
      {
         _ar & itermax;
         _ar & tolerance;
         _ar & linetolerance;
         _ar & linestep;
      }
    
    //! Interface for evaluating the functioanl
    template<typename T_FUNTIONAL>
    double gsl_func_f( const gsl_vector *_x, void *_params)
    {
      gsl_const_iterator i_x(_x);
      gsl_const_iterator i_x_end(i_x); i_x_end += ((T_FUNTIONAL*) _params)->size();
      typename T_FUNTIONAL::iterator i_var = ((T_FUNTIONAL*) _params)->begin();
      for(; i_x != i_x_end; ++i_x, ++i_var)
        *i_var = *i_x;

      return ((T_FUNTIONAL*) _params)->evaluate();
    }
    //! Interface for evaluating the derivatives of functional
    template<typename T_FUNTIONAL>
    void gsl_func_df( const gsl_vector *_x, void *_params, gsl_vector *_grad)
    {
      gsl_const_iterator i_x(_x);
      gsl_const_iterator i_x_end(i_x); i_x_end += ((T_FUNTIONAL*) _params)->size();
      typename T_FUNTIONAL::iterator i_var = ((T_FUNTIONAL*) _params)->begin();
      for(; i_x != i_x_end; ++i_x, ++i_var)
        *i_var = *i_x;
      ((T_FUNTIONAL*) _params)->evaluate_gradient(_grad->data);
    }
    //! Interface for evaluating the functional and its derivatives.
    template<typename T_FUNTIONAL>
    void gsl_func_fdf( const gsl_vector *_x, void *_params, double *_r, gsl_vector *_grad)
    {
      gsl_const_iterator i_x(_x);
      gsl_const_iterator i_x_end(i_x); i_x_end += ((T_FUNTIONAL*) _params)->size();
      typename T_FUNTIONAL::iterator i_var = ((T_FUNTIONAL*) _params)->begin();
      for(; i_x != i_x_end; ++i_x, ++i_var)
        *i_var = *i_x;

      *_r = ((T_FUNTIONAL*) _params)->evaluate_with_gradient(_grad->data);
    }
  }
} // namespace LaDa
#endif
