//
//  Version: $Id$
//
#ifndef _OPT_MINIMIZE_GSL_H_
#define _OPT_MINIMIZE_GSL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/opt_minimize_base.h>
#include <gsl/gsl_multimin.h>
#include <functional>
#include <algorithm>
#include <opt/compose_functors.h>

#ifdef _MPI
#  include "mpi/mpi_object.h"
#endif

namespace minimizer {

  class gsl_iterator
  {
    public:
      typedef double* iterator_category;
      typedef double value_type;
      typedef int difference_type;
      typedef double* pointer;
      typedef double reference;

    private:
      gsl_vector *x;
      double* iterator;

    public:
      gsl_iterator( gsl_vector *_x ) : x(_x), iterator(x->data) {};
      gsl_iterator( const gsl_iterator &_x ) : x(_x.x), iterator(_x.iterator) {};
      gsl_iterator& operator++()
       { iterator += x->stride; return *this; }
      gsl_iterator& operator++(int)
       { iterator += x->stride; return *this; }
      gsl_iterator& operator+=(types::t_unsigned n)
       { iterator += n*x->stride; return *this; }
      gsl_iterator operator+(types::t_unsigned n)
       { gsl_iterator it(x); it += n; return it;}
      gsl_iterator operator-(types::t_unsigned n)
       { gsl_iterator it(x); it -= n; return it;}
      gsl_iterator& operator--()
       { iterator -= x->stride; return *this; }
      gsl_iterator& operator-=(types::t_unsigned n)
       { iterator -= n*x->stride; return *this; }
      double& operator*() const
       { return *iterator; }
      bool operator==(const gsl_iterator &_x) const
       { return iterator == _x.iterator; }
      bool operator!=(const gsl_iterator &_x) const
       { return iterator != _x.iterator; }
  };
  class gsl_const_iterator
  {
    public:
      typedef const double* iterator_category;
      typedef const double value_type;
      typedef int difference_type;
      typedef const double* pointer;
      typedef const double reference;

    private:
      const gsl_vector *x;
      pointer iterator;
    public:
      gsl_const_iterator( const gsl_vector *_x ) : x(_x), iterator(x->data) {};
      gsl_const_iterator( const gsl_const_iterator &_x ) : x(_x.x), iterator(_x.iterator) {};
      gsl_const_iterator& operator++()
       { iterator += x->stride; return *this; }
      gsl_const_iterator& operator++(int)
       { iterator += x->stride; return *this; }
      gsl_const_iterator& operator+=(types::t_unsigned n)
       { iterator += n*x->stride; return *this; }
      gsl_const_iterator operator+(types::t_unsigned n)
       { gsl_const_iterator it(x); it += n; return it;}
      gsl_const_iterator operator-(types::t_unsigned n)
       { gsl_const_iterator it(x); it -= n; return it;}
      gsl_const_iterator& operator--()
       { iterator -= x->stride; return *this; }
      gsl_const_iterator& operator-=(types::t_unsigned n)
       { iterator -= n*x->stride; return *this; }
      reference& operator*() const
       { return *iterator; }
      bool operator==(const gsl_const_iterator &_x)
       { return iterator == _x.iterator; }
      bool operator!=(const gsl_const_iterator &_x)
       { return iterator != _x.iterator; }
  };

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
    
  template<typename T_FUNCTIONAL> 
  class GnuSL : public Base< T_FUNCTIONAL >
  {
    public:
      typedef T_FUNCTIONAL t_Functional;
    protected:
      using Base<t_Functional>::current_func;
      types::t_real tolerance;
      types::t_real linetolerance;
      types::t_real linestep;
      types::t_unsigned iter; // number of performed iterations
      types::t_unsigned itermax; // maximum number of iterations
      enum t_gsl_minimizer_type { GSL_NONE, GSL_FR, GSL_PR,                       
                                  GSL_BFGS2, GSL_BFGS, GSL_SD };
      t_gsl_minimizer_type type; // maximum number of iterations
    public:
      bool do_print;
                                  
      
    public:
      GnuSL () : tolerance(types::tolerance),
                 linetolerance(0.01),
                 linestep(0.1), iter(0),
                 itermax(500), do_print(false) {}
      explicit
        inline GnuSL   (t_Functional& _func)
                     : Base<t_Functional>( _func ),
                       tolerance(types::tolerance),
                       linetolerance(0.01),
                       linestep(0.1), iter(0),
                       itermax(500), do_print(false) {}
      virtual ~GnuSL(){};

      virtual bool minimize()
      {
        // initializes object related stuff
        if( not current_func->init() )
          return false;
        int status;
        
        const gsl_multimin_fdfminimizer_type *T;
        gsl_multimin_fdfminimizer *s;
        
        gsl_vector *x;
        
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
          case GSL_FR: T = gsl_multimin_fdfminimizer_conjugate_fr; break;
          case GSL_PR: T = gsl_multimin_fdfminimizer_conjugate_pr; break;
          case GSL_BFGS: T = gsl_multimin_fdfminimizer_vector_bfgs; break;
          case GSL_BFGS2: T = gsl_multimin_fdfminimizer_vector_bfgs2; break;
          case GSL_SD: T = gsl_multimin_fdfminimizer_steepest_descent; break;
        }
        s = gsl_multimin_fdfminimizer_alloc (T, gsl_func.n);
        if (not s)
          return false;
        
        gsl_multimin_fdfminimizer_set (s, &gsl_func, x, linestep, linetolerance);
        
        double newe=current_func->energy(), olde=0;
        types::t_unsigned nomove=0;
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
          if ( tolerance > std::abs(newe - olde) )
            ++nomove;
          else 
            nomove = 0;
          if ( do_print )
            std::cout << current_func->energy() << " " << iter << std::endl;
        
        }
        while ( nomove < 3 and iter < itermax);
        
        
        gsl_func_f<t_Functional>( gsl_multimin_fdfminimizer_x(s), (void*)current_func );
        if ( do_print )
          std::cout << current_func->energy() << " " << iter << " " << nomove << std::endl;

        gsl_multimin_fdfminimizer_free (s);
        gsl_vector_free (x);

        return status == GSL_SUCCESS;
      }  // dummy minimizer


      bool Load( const TiXmlElement &_node )
      {
        const TiXmlElement *parent;
        std::string str;
      
        // This whole section tries to find a <Functional type="vff"> tag
        // in _node or its child
        str = _node.Value();
        if ( str.compare("Minimizer" ) != 0 )
          parent = _node.FirstChildElement("Minimizer");
        else
          parent = &_node;
      
        
        while (parent)
        {
          str = "";
          type = GSL_NONE;
          if ( parent->Attribute( "type" )  )
            str = parent->Attribute("type");
          if ( str.compare("gsl_fr" ) == 0 )
            type = GSL_FR;
          if ( str.compare("gsl_pr" ) == 0 )
            type = GSL_PR;
          if ( str.compare("gsl_bfgs2" ) == 0 )
            type = GSL_BFGS2;
          if ( str.compare("gsl_bfgs" ) == 0 )
            type = GSL_BFGS;
          if ( str.compare("gsl_sd" ) == 0 )
            type = GSL_SD;
          if ( type != GSL_NONE )
            break;
          parent = parent->NextSiblingElement("Minimizer");
        }
        if ( not parent )
        {
          std::cerr << "Could not find an <Minimizer type=\"some gsl\"> tag in input file" 
                    << std::endl;
          return false;
        } 
      
        double d; int n;
        parent->Attribute( "itermax", &n );
        itermax = (n > 0) ? abs(n) : 10000;
        parent->Attribute( "tolerance", &d );
        tolerance = d ? types::t_real(d) : types::tolerance;
        parent->Attribute( "linetolerance", &d );
        linetolerance = d ? types::t_real(d) : 0.01;
        parent->Attribute( "linestep", &d );
        linestep = d ? types::t_real(d) : 0.1;
      
        return true;
     }

#ifdef _MPI
     bool broadcast( mpi::BroadCast &_bc )
     {
        if ( not _bc.serialize( itermax ) ) return false;
        if ( not _bc.serialize( tolerance ) ) return false;
        if ( not _bc.serialize( linetolerance ) ) return false;
        return _bc.serialize( linestep );
     }
#endif 
  };

}

#endif
