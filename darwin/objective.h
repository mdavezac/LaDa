//
//  Version: $Id$
//
#ifndef _MULTIOB_OBJECTIVE_H_
#define _MULTIOB_OBJECTIVE_H_

#include<vector>
#include<math.h>
#include<stdexcept>

#include <tinyxml/tinyxml.h>

#include "opt/types.h"
#include "opt/convex_hull.h"
#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

#include "gatraits.h"
#include "loadsave.h"
#include "print_xmgrace.h"

namespace Objective
{
  template< class T_EVALUATOR, class T_GA_TRAITS = Traits::GA<T_EVALUATOR>, 
            class T_QUANTITY_TRAITS = typename T_GA_TRAITS :: t_QuantityTraits,
            class T_VA_TRAITS = typename T_GA_TRAITS :: t_VA_Traits >
  class Base
  {
    public: 
      typedef T_EVALUATOR t_Evaluator;
      typedef T_GA_TRAITS t_GA_Traits;
      typedef T_QUANTITY_TRAITS t_QuantityTraits;
      typedef T_VA_TRAITS t_VA_Traits;
    protected:
      typedef typename t_GA_Traits :: t_Individual t_Individual;
      typedef typename t_QuantityTraits :: t_Quantity t_Quantity;
      typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type t_VA_Type;
      typedef darwin::SaveObject<t_Evaluator> t_SaveOp;
      typedef darwin::LoadObject<t_Evaluator> t_LoadOp;
    protected:
      static const t_Individual *current_indiv;
    public:
      bool const static is_scalar = t_QuantityTraits :: is_scalar;
      bool const static is_vector = t_QuantityTraits :: is_vectorial;
    public:
      Base() {}
      virtual ~Base() {}

      void static init(const t_Individual& _indiv);
      virtual t_ScalarQuantity operator()( const t_Quantity& ) = 0;
      virtual void evaluate_gradient( const t_Quantity &,
                                      t_QuantityGradients &,
                                      t_ScalarQuantity *) = 0;
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &,
                                                       t_QuantityGradients&,
                                                       t_VA_Type *) = 0;
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                               t_QuantityGradients&,
                                               types::t_unsigned) = 0;
      virtual bool is_valid() const = 0;
      virtual std::string what_is() const = 0;
      virtual std::string print() const = 0;
      virtual bool Save( TiXmlElement &_node, t_SaveOp& _op) { return true; };
      virtual bool Restart( const  TiXmlElement &_node, t_LoadOp &_op) { return true; };
  };
  template< class T_EVALUATOR, class T_GA_TRAITS, class T_QUANTITY_TRAITS, class T_VA_TRAITS >
  void Base<T_EVALUATOR, T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: init( const t_Individual & _indiv)
    { Base<T_EVALUATOR, T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: current_indiv = &_indiv; }
  template< class T_EVALUATOR, class T_GA_TRAITS, class T_QUANTITY_TRAITS, class T_VA_TRAITS >
    const typename Base<T_EVALUATOR, T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: t_Individual* 
      Base<T_EVALUATOR, T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: current_indiv = NULL;
  
  template < class T_EVALUATOR, class T_GA_TRAITS = Traits::GA<T_EVALUATOR> >
    struct Types
    {
      public:
        typedef T_EVALUATOR t_Evaluator;
        typedef T_GA_TRAITS t_GA_Traits;
      protected:
        typedef typename t_Evaluator :: t_Individual t_Individual;
        typedef typename t_Evaluator :: t_IndivTraits t_IndivTraits;
        typedef typename t_IndivTraits :: t_QuantityTraits t_QuantityTraits;
        typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;
      public:
      typedef Base< t_Evaluator, t_GA_Traits, 
                    typename Traits::Quantity< t_ScalarQuantity > > Scalar;
      typedef Base<t_Evaluator, t_GA_Traits >  Vector;
      static Vector* new_from_xml( const TiXmlElement &_node );
    };

  template< class T_EVALUATOR, class T_GA_TRAITS >
  class Maximize : public Types< T_EVALUATOR, T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_EVALUATOR t_Evaluator;
      typedef T_GA_TRAITS t_GA_Traits;
    protected:
      typedef typename Types<t_Evaluator, t_GA_Traits> :: Scalar t_Base;
      typedef typename t_GA_Traits :: t_Individual t_Individual;
      typedef typename t_GA_Traits :: t_IndivTraits t_IndivTraits;
      typedef typename t_Base :: t_Quantity t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity t_ScalarQuantity;
      typedef typename t_Base :: t_VA_Traits t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type t_VA_Type;

    public:
      Maximize() {}
      Maximize( const Maximize &) {}
      virtual ~Maximize() {}
      
      virtual t_ScalarQuantity operator()(const t_Quantity& _val) { return -_val; }
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad)
      { 
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
          *i_grad_result -= *i_grad;
      }
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad)  
      { 
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
          *i_grad_result -= *i_grad;
        return -_val;
      }
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n) { return -_grad[_n]; }
      bool is_valid() const { return true; }
      virtual std::string what_is() const { return " Maximize"; }
      virtual std::string print() const { return ""; }
  };
  template< class T_EVALUATOR, class T_GA_TRAITS >
  class Minimize : public Types< T_EVALUATOR, T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_EVALUATOR t_Evaluator;
      typedef T_GA_TRAITS t_GA_Traits;
    protected:
      typedef typename Types<t_Evaluator, t_GA_Traits> :: Scalar t_Base;
      typedef typename t_GA_Traits :: t_Individual t_Individual;
      typedef typename t_GA_Traits :: t_IndivTraits t_IndivTraits;
      typedef typename t_Base :: t_Quantity t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity t_ScalarQuantity;
      typedef typename t_Base :: t_VA_Traits t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type t_VA_Type;

    public:
      Minimize() {}
      Minimize( const Minimize &) {}
      virtual ~Minimize() {}
      
      virtual t_ScalarQuantity operator()(const t_Quantity& _val) { return _val; }
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad)
      { 
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
          *i_grad_result += *i_grad;
      }
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad)  
      { 
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
          *i_grad_result += *i_grad;
        return _val;
      }
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n) { return _grad[_n]; }
      bool is_valid() const { return true; }
      virtual std::string what_is() const { return " Minimize"; }
      virtual std::string print() const { return ""; }
  };
  template< class T_EVALUATOR, class T_GA_TRAITS >
  class Target : public Types< T_EVALUATOR, T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_EVALUATOR t_Evaluator;
      typedef T_GA_TRAITS t_GA_Traits;
    protected:
      typedef typename Types<t_Evaluator, t_GA_Traits> :: Scalar t_Base;
      typedef typename t_GA_Traits :: t_Individual t_Individual;
      typedef typename t_GA_Traits :: t_IndivTraits t_IndivTraits;
      typedef typename t_Base :: t_Quantity t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity t_ScalarQuantity;
      typedef typename t_Base :: t_VA_Traits t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type t_VA_Type;

    protected:
      t_ScalarQuantity target; 
    public:
      Target( t_ScalarQuantity _target ) : target( _target ) {}
      Target( const Target &_c ) : target( _c.target ) {}
      virtual ~Target() {}
      
      virtual t_ScalarQuantity operator()(const t_Quantity& _val)
        { return std::abs( _val - target ); }
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad)
      {
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        if ( _val > target ) 
          for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
            *i_grad_result += *i_grad;
        else 
          for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
            *i_grad_result -= *i_grad;
      }
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad)  
      {
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        if ( _val > target ) 
          for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
            *i_grad_result += *i_grad;
        else 
          for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
            *i_grad_result -= *i_grad;
        return std::abs(target - _val);
      }
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity & _val,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n) 
      {
        return ( _val > target ) ? _grad[_n]: -_grad[_n];
      }
      bool is_valid() const { return true; }
      virtual std::string print() const { return ""; }
      virtual std::string what_is() const
      { 
        std::ostringstream sstr;
        sstr << " Target (" << target << ")";
        return sstr.str();
      }
  };
  template< class T_EVALUATOR, class T_GA_TRAITS >
  class ConvexHull : public Types< T_EVALUATOR, T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_EVALUATOR t_Evaluator;
      typedef T_GA_TRAITS t_GA_Traits;
    protected:
      typedef typename Types<t_Evaluator, t_GA_Traits> :: Scalar t_Base;
      typedef typename t_GA_Traits :: t_Individual t_Individual;
      typedef typename t_GA_Traits :: t_IndivTraits t_IndivTraits;
      typedef typename t_Base :: t_Quantity t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity t_ScalarQuantity;
      typedef opt::ConvexHull::Base<t_Individual> t_ConvexHull;
      typedef typename t_Base :: t_VA_Traits t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type t_VA_Type;
      typedef darwin::SaveObject<t_Evaluator> t_SaveOp;
      typedef darwin::LoadObject<t_Evaluator> t_LoadOp;

    protected:
      t_ConvexHull convexhull; 
      mutable bool valid;
      using t_Base :: current_indiv;

    public:
      ConvexHull() : valid(true) {}
      ConvexHull( const ConvexHull &_c ) : convexhull( _c.convexhull ), valid(true) {}
      virtual ~ConvexHull() {}
      
      virtual t_ScalarQuantity operator()(const t_Quantity& _val);
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad)
        { evaluate_with_gradient( _val, _grad, _i_grad ); }

      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *);
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity & _val,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n);
      bool is_valid() const
      {
        if ( valid ) return true;
        valid = true;
        return false;
      }
      virtual std::string what_is() const { return " Convex-Hull"; }
      virtual std::string print() const { return convexhull.print(); }
      virtual bool Save( TiXmlElement &_node, t_SaveOp& _op)
        { return convexhull.Save( _node, _op ); };
      virtual bool Restart( const  TiXmlElement &_node, t_LoadOp &_op)
        { return convexhull.Load( _node, _op ); };
  };
  template< class T_EVALUATOR, class T_GA_TRAITS >
  typename ConvexHull<T_EVALUATOR, T_GA_TRAITS> :: t_ScalarQuantity
    ConvexHull<T_EVALUATOR, T_GA_TRAITS> :: operator()(const t_Quantity& _val)
    {
      t_Quantity x = current_indiv->get_concentration();
      t_Quantity base = (t_Quantity) convexhull.evaluate( x );
    
      if ( _val >= base ) return _val - base;
    
      if ( convexhull.add( _val, *current_indiv ) )
        valid = false;

      return 0.0;
    }
  template< class T_EVALUATOR, class T_GA_TRAITS >
  typename ConvexHull<T_EVALUATOR, T_GA_TRAITS> :: t_ScalarQuantity
    ConvexHull<T_EVALUATOR, T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_val,
                                                                        t_QuantityGradients &_grad,
                                                                        t_VA_Type *_i_grad)
    {
      t_Quantity x = current_indiv->get_concentration();
      t_Quantity base = (t_Quantity) convexhull.evaluate( x );
      types::t_real Ninv = 1.0 / ( (types::t_real ) current_indiv->Object().Container().size() );
      types::t_real gradient = convexhull.evaluate_gradient( x ) * Ninv;
      typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
      typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
      t_VA_Type *i_grad_result = _i_grad;
      for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
        *i_grad_result += *i_grad - gradient;
      
      if ( _val >= base ) return _val - base;

      if ( convexhull.add( _val, *current_indiv ) )
        valid = false;

      return 0.0;
    }
  template< class T_EVALUATOR, class T_GA_TRAITS >
  typename ConvexHull<T_EVALUATOR, T_GA_TRAITS> :: t_VA_Type
    ConvexHull<T_EVALUATOR, T_GA_TRAITS> :: evaluate_one_gradient( const t_Quantity &_val,
                                                                   t_QuantityGradients &_grad,
                                                                   types::t_unsigned _n) 
    {
      t_Quantity x = current_indiv->get_concentration();
      types::t_real Ninv = 1.0 / ( (types::t_real ) current_indiv->Object().Container().size() );
      types::t_real gradient = convexhull.evaluate_gradient( x ) * Ninv;
      return _grad[_n] - gradient;
    }


  template< class T_EVALUATOR, class T_GA_TRAITS >
  class Container : public Types<T_EVALUATOR, T_GA_TRAITS> :: Vector
  {
    public:
      typedef T_EVALUATOR t_Evaluator;
      typedef T_GA_TRAITS t_GA_Traits;
    protected:
      typedef typename Types<t_Evaluator, t_GA_Traits> :: Vector t_Base;
      typedef typename t_GA_Traits :: t_Individual t_Individual;
      typedef typename t_GA_Traits :: t_IndivTraits t_IndivTraits;
      typedef typename t_Base :: t_Quantity t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity t_ScalarQuantity;
      typedef typename Types<t_Individual, t_IndivTraits> :: Scalar t_Objective;
      typedef std::vector< t_Objective > t_Objectives;
      typedef typename t_Base :: t_VA_Traits t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type t_VA_Type;

    protected:
      t_Objectives objectives;

    public:
      Container() {}
      Container( const Container &_c ) : objectives(_c.objectives) {}
      virtual ~Container () 
      {
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.end();
        for(; i_objective != i_end; ++i_objective )
          delete *i_objective;
        objectives.clear();
      }

      bool is_valid() const
      {
        typename t_Objectives::const_iterator i_objective = objectives.begin();
        typename t_Objectives::const_iterator i_end = objectives.begin();
        for(; i_objective != i_end; ++i_objective )
          if ( not (*i_objective)->is_valid() )
            return false;
        return true;
      }
  };

  template< class T_EVALUATOR, class T_GA_TRAITS >
  class LinearSum : public Container<T_EVALUATOR, T_GA_TRAITS>
  {
    public:
      typedef T_EVALUATOR t_Evaluator;
      typedef T_GA_TRAITS t_GA_Traits;
    protected:
      typedef Container<t_Evaluator, t_GA_Traits> t_Base;
      typedef typename t_GA_Traits :: t_Individual t_Individual;
      typedef typename t_GA_Traits :: t_IndivTraits t_IndivTraits;
      typedef typename t_IndivTraits :: t_QuantityTraits t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_Quantity t_Quantity;
      typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;
      typedef typename Types<t_Individual, t_IndivTraits> :: Scalar t_Objective;
      typedef std::vector< t_Objective > t_Objectives;

    protected:
      using t_Base :: objectives;
      std::vector< t_ScalarQuantity > coefs;

    public:
      LinearSum() {}
      LinearSum ( const LinearSum &_c ) : t_Base(_c), coefs(_c.coefs) {}
      virtual ~LinearSum() {}

      virtual t_ScalarQuantity operator()( const t_Quantity& _val ) 
      {
        if ( _val.size() != coefs.size() )
          throw std::runtime_error( "Wrong number of objective functions\n" );

        types::t_real inter = 0;
        typename t_Quantity :: iterator i_val = _val.begin();
        typename t_Quantity :: iterator i_val_end = _val.end();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.begin();
        for(; i_objective != i_end and i_val != i_val_end; ++i_objective, ++i_coef, ++i_val )
          inter += ( *i_coef ) * (*i_objective)()( *i_val );
          
        return inter;
      };
  };

  template< class T_EVALUATOR, class T_GA_TRAITS >
  typename Types<T_EVALUATOR, T_GA_TRAITS> :: Vector* 
     Types<T_EVALUATOR, T_GA_TRAITS> :: new_from_xml( const TiXmlElement &_node )
      {
        if ( not &_node ) return NULL;
        std::string str = "minimize"; 
        std::string name = _node.Value();
        if (    name.compare("Objective") == 0 
             or name.compare("Method") == 0 )
        {
          if ( _node.Attribute( "type" ) )
            str = _node.Attribute( "type" );
        }
        else if ( _node.Attribute("objective") )
          str = _node.Attribute( "objective" );
        if ( str.compare("convexhull") == 0 )
        {
          darwin::printxmg << darwin::PrintXmg::comment << "Objective: ConvexHull" << darwin::PrintXmg::endl;
          return new ConvexHull<t_Evaluator, t_GA_Traits>;
        }
        else if ( str.compare("minimize") == 0 )
        {
          darwin::printxmg << darwin::PrintXmg::comment << "Objective: Minimize" << darwin::PrintXmg::endl;
          return new Minimize<t_Evaluator, t_GA_Traits>;
        }
        else if ( str.compare("maximize") == 0 )
        {
          darwin::printxmg << darwin::PrintXmg::comment << "Objective: Maximize" << darwin::PrintXmg::endl;
          return new Maximize<t_Evaluator, t_GA_Traits>;
        }
        else if (str.compare("target") == 0 )
        {
          if( _node.Attribute("target") )
          {
            double d; _node.Attribute("target", &d );
            darwin::printxmg << darwin::PrintXmg::comment
                             << "Objective: Target (" << d << ")" << darwin::PrintXmg::endl;
            return new Target<t_Evaluator, t_GA_Traits>( (types::t_real) d );
          }
        }
        if ( _node.FirstChildElement( "Objective" ) )
         return new_from_xml( *_node.FirstChildElement( "Objective" ) ); 

        return NULL;
      }

}  // namespace fitness

namespace darwin
{
  class Fitness
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<darwin::Fitness>( darwin::Fitness & ); 
#endif 
    public:
      typedef types::t_real t_Quantity;

    protected:
      t_Quantity quantity;
      bool is_valid;

    public:
      Fitness() : is_valid( false )  {}
      Fitness( const Fitness & _c ) : quantity( _c.quantity ), is_valid( _c.is_valid ) {}
      ~Fitness() {}


      bool operator<(const Fitness & _f) const
        { return std::abs(quantity - _f.quantity) > types::tolerance and quantity > _f.quantity; }
      bool operator>(const Fitness & _f) const
        { return std::abs(quantity - _f.quantity) > types::tolerance and quantity < _f.quantity; }
      bool operator==(const Fitness & _f) const
        { return std::abs(quantity - _f.quantity) < types::tolerance; }

      bool invalid() const { return not is_valid; }
      void invalidate() { is_valid = false; }
      void operator=( types::t_real _fit ) 
      {
        is_valid = true;
        quantity = _fit;
      }
      operator t_Quantity() const
      { 
        if ( not is_valid )
          throw std::runtime_error( " Invalid Fitness !!\n" );
        return quantity;
      }
      bool Load( const TiXmlElement & _node );
      bool Save( TiXmlElement & _node ) const;
  };

}

std::ostream & operator<<( std::ostream &_os, darwin::Fitness &_fit );

#endif //  _MULTIOB_OBJECTIVE_H_
