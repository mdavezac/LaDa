//
//  Version: $Id$
//
#ifndef _OBJECTIVE_IMPL_H_
#define _OBJECTIVE_IMPL_H_

#include <opt/debug.h>

namespace GA
{

namespace Objective
{
  template< class T_GA_TRAITS, class T_QUANTITY_TRAITS, class T_VA_TRAITS >
  void Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: init( const t_Individual & _indiv)
    { Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: current_indiv = &_indiv; }
  template<class T_GA_TRAITS, class T_QUANTITY_TRAITS, class T_VA_TRAITS >
    const typename Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: t_Individual* 
      Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: current_indiv = NULL;
  template<class T_GA_TRAITS, class T_QUANTITY_TRAITS, class T_VA_TRAITS >
    typename Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: t_Fitness 
      Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: fitness;

  template< class T_GA_TRAITS >
  inline void Maximize<T_GA_TRAITS> :: evaluate_gradient( const t_Quantity &_val,
                                                          t_QuantityGradients &_grad,
                                                          t_VA_Type *_i_grad)
  { 
    typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
    typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
    t_VA_Type *i_grad_result = _i_grad;
    for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
      *i_grad_result -= *i_grad;
  }
  template< class T_GA_TRAITS >
  inline typename Maximize<T_GA_TRAITS>::t_ScalarQuantity
    Maximize<T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_val,
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

  template< class T_GA_TRAITS >
  inline void Minimize<T_GA_TRAITS> :: evaluate_gradient( const t_Quantity &_val,
                                                          t_QuantityGradients &_grad,
                                                          t_VA_Type *_i_grad)
  { 
    typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
    typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
    t_VA_Type *i_grad_result = _i_grad;
    for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
      *i_grad_result += *i_grad;
  }

  template< class T_GA_TRAITS >
  inline typename Minimize<T_GA_TRAITS>::t_ScalarQuantity
    Minimize<T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_val,
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


  template< class T_GA_TRAITS >
  inline void Target< T_GA_TRAITS > :: evaluate_gradient( const t_Quantity &_q,
                                                          t_QuantityGradients &_grad,
                                                          t_VA_Type *_i_grad)
  {
    typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
    typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
    t_VA_Type *i_grad_result = _i_grad;
    if ( _q > q_0 ) 
      for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
        *i_grad_result += *i_grad;
    else 
      for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
        *i_grad_result -= *i_grad;
  }
  template< class T_GA_TRAITS >
  inline typename Target<T_GA_TRAITS> :: t_ScalarQuantity
      Target< T_GA_TRAITS > :: evaluate_with_gradient( const t_Quantity &_q,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad)  
  {
    typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
    typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
    t_VA_Type *i_grad_result = _i_grad;

    if ( t_QuantityTraits::gt(_q, q_0 ) )
      for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
        *i_grad_result += *i_grad;
    else if ( t_QuantityTraits::le(_q, q_0 ) )
      for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
        *i_grad_result -= *i_grad;
    return std::abs(q_0 - _q);
  }
  template< class T_GA_TRAITS >
  inline typename Target<T_GA_TRAITS> :: t_VA_Type
      Target< T_GA_TRAITS > :: evaluate_one_gradient( const t_Quantity & _q,
                                                      t_QuantityGradients& _grad,
                                                      types::t_unsigned _n) 
  {
    if ( t_QuantityTraits::gt(_q, q_0 ) ) return _grad[_n];
    else if ( t_QuantityTraits::le(_q, q_0 ) ) return -_grad[_n];
    return t_VA_Type(0);
  }
  template< class T_GA_TRAITS >
  inline std::string Target< T_GA_TRAITS > :: what_is() const
  { 
    std::ostringstream sstr;
    sstr << " Target (" << q_0 << ")";
    return sstr.str();
  }





  template< class T_GA_TRAITS >
  inline bool ConvexHull<T_GA_TRAITS> :: is_valid() const
  {
    if ( valid ) return true;
    valid = true;
    return false;
  }
  template< class T_GA_TRAITS >
  const typename ConvexHull<T_GA_TRAITS> :: t_Fitness&
    ConvexHull<T_GA_TRAITS> :: operator()(const t_Quantity& _q)
    {
      t_Quantity x = current_indiv->get_concentration();
      t_Quantity base = (t_Quantity) convexhull.evaluate( x );
    
      if ( _q >= base and convexhull.size() >= 2 )
       { fitness = _q - base;  return fitness; }
    
      if ( convexhull.add( _q, *current_indiv ) )
        valid = false;

      fitness = t_ScalarQuantity(0);
      return fitness;
    }
  template< class T_GA_TRAITS >
  void ConvexHull<T_GA_TRAITS> :: evaluate_gradient( const t_Quantity &_q,
                                                     t_QuantityGradients &_grad,
                                                     t_VA_Type *_i_grad)
  {
    t_Quantity x = current_indiv->get_concentration();

    typedef typename t_GATraits :: t_Object :: t_Container t_Container;
    const t_Container &container = current_indiv->Object().Container(); 
    types::t_real Ninv = 1.0 / ( (types::t_real ) container.size() );
    types::t_real left = convexhull.evaluate_left_gradient( x ) * Ninv;
    types::t_real right = convexhull.evaluate_right_gradient( x ) * Ninv;
    typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
    typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
    typename t_Container :: const_iterator i_var = container.begin(); 
    t_VA_Type *i_grad_result = _i_grad;
    for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result, ++i_var )
      *i_grad_result += *i_grad - (*i_var > 0e0 ? left: right);
  }
  template< class T_GA_TRAITS >
  typename ConvexHull<T_GA_TRAITS> :: t_ScalarQuantity
    ConvexHull<T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_q,
                                                       t_QuantityGradients &_grad,
                                                       t_VA_Type *_i_grad)
    {
      t_Quantity x = current_indiv->get_concentration();
      t_Quantity base = (t_Quantity) convexhull.evaluate( x );

      typedef typename t_GATraits :: t_Object :: t_Container t_Container;
      const t_Container &container = current_indiv->Object().Container(); 
      types::t_real Ninv = 1.0 / ( (types::t_real ) container.size() );
      types::t_real left = convexhull.evaluate_left_gradient( x ) * Ninv;
      types::t_real right = convexhull.evaluate_right_gradient( x ) * Ninv;
      typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
      typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
      typename t_Container :: const_iterator i_var = container.begin(); 
      t_VA_Type *i_grad_result = _i_grad;
      for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result, ++i_var )
        *i_grad_result += *i_grad - (*i_var > 0e0 ? left: right);
      
      if ( _q >= base ) return _q - base;

      if ( convexhull.add( _q, *current_indiv ) )
        valid = false;

      return t_ScalarQuantity(0);
    }
  template< class T_GA_TRAITS >
  typename ConvexHull<T_GA_TRAITS> :: t_VA_Type
    ConvexHull<T_GA_TRAITS> :: evaluate_one_gradient( const t_Quantity &_q,
                                                      t_QuantityGradients &_grad,
                                                      types::t_unsigned _n) 
    {
      t_Quantity x = current_indiv->get_concentration();
      types::t_real Ninv = 1.0 / ( (types::t_real ) current_indiv->Object().Container().size() );
      types::t_real gradient 
         = Ninv * ( ( current_indiv->Object().Container()[_n] > 0e0 ) ? 
                       convexhull.evaluate_left_gradient( x ):
                       convexhull.evaluate_right_gradient( x ) );
      return _grad[_n] - gradient;
    }






  template<class T_GA_TRAITS >
  Container<T_GA_TRAITS> :: ~Container () 
  {
    typename t_Objectives :: iterator i_objective = objectives.begin();
    typename t_Objectives :: iterator i_end = objectives.end();
    for(; i_objective != i_end; ++i_objective ) delete *i_objective;
    objectives.clear();
  }

  template<class T_GA_TRAITS >
  inline bool Container<T_GA_TRAITS> :: is_valid() const
  {
    typename t_Objectives :: const_iterator i_objective = objectives.begin();
    typename t_Objectives :: const_iterator i_end = objectives.end();
    for(; i_objective != i_end; ++i_objective )
      if ( not (*i_objective)->is_valid() ) return false;
    return true;
  }
  template<class T_GA_TRAITS >
  inline bool Container<T_GA_TRAITS> :: Save( TiXmlElement &_node, t_SaveOp& _op)
  {
    typename t_Objectives :: const_iterator i_objective = objectives.begin();
    typename t_Objectives :: const_iterator i_end = objectives.end();
    for(; i_objective != i_end; ++i_objective )
      if ( not (*i_objective)->Save( _node, _op ) ) return false;
    return true;
  }
  template<class T_GA_TRAITS >
  inline bool Container<T_GA_TRAITS> :: Restart( const  TiXmlElement &_node, t_LoadOp &_op) 
  {
    typename t_Objectives :: iterator i_objective = objectives.begin();
    typename t_Objectives :: iterator i_end = objectives.end();
    for(; i_objective != i_end; ++i_objective )
      if ( not (*i_objective)->Restart( _node, _op ) ) return false;
    return true;
  }
  template<class T_GA_TRAITS >
  inline bool Container<T_GA_TRAITS> :: does_store() const
  {
    typename t_Objectives :: const_iterator i_objective = objectives.begin();
    typename t_Objectives :: const_iterator i_end = objectives.end();
    for(; i_objective != i_end; ++i_objective )
      if ( (*i_objective)->does_store() ) return true;
    return false;
  }
  template<class T_GA_TRAITS >
  inline std::string Container<T_GA_TRAITS> :: print() const
  {
    typename t_Objectives :: const_iterator i_objective = objectives.begin();
    typename t_Objectives :: const_iterator i_end = objectives.end();
    std::ostringstream sstr;
    for(; i_objective != i_end; ++i_objective )
      sstr << (*i_objective)->print();
    return sstr.str();
  }
















  template<class T_GA_TRAITS >
  inline void LinearSum<T_GA_TRAITS> :: add( t_Objective *_objective, t_ScalarQuantity _coef )
  {
    if ( not _objective ) return;
    coefs.push_back( _coef );
    objectives.push_back( _objective );
  }
  template<class T_GA_TRAITS >
  inline std::string LinearSum<T_GA_TRAITS> :: what_is() const
  {
    std::ostringstream sstr;
    sstr << "LinearSum begin{ ";
    typename t_Objectives::const_iterator i_objective = objectives.begin();
    typename t_Objectives::const_iterator i_end = objectives.end();
    typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
    for(; i_objective != i_end; ++i_objective, ++i_coef )
      sstr << (*i_objective)->what_is() << "[" << *i_coef << "] ";
    sstr << " }end"; 
    return  sstr.str();
  }
  template< class T_GA_TRAITS >
    const typename LinearSum<T_GA_TRAITS>::t_Fitness&
      LinearSum<T_GA_TRAITS> :: operator()( const t_Quantity& _q ) 
      {
        __DOASSERT( t_QuantityTraits::size(_q) != coefs.size(),
                       "Wrong number of objective functions during evaluation \n"
                    << "Expected " << t_QuantityTraits::size(_q) << "\n"
                    << "Found " << coefs.size() << "\n" )

        t_ScalarQuantity inter = 0;
        typename t_Quantity :: const_iterator i_val = _q.begin();
        typename t_Quantity :: const_iterator i_val_end = _q.end();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.end();
        fitness.clear();
        for(; i_objective != i_end and i_val != i_val_end; ++i_objective, ++i_coef, ++i_val )
        {
          double r = (*i_objective)->operator()( *i_val );
          fitness.push_back( r );
          inter += ( *i_coef ) * r;
        }
          
        fitness = inter;
        return fitness;
      };
  template< class T_GA_TRAITS >
    typename LinearSum<T_GA_TRAITS>::t_ScalarQuantity
      LinearSum<T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_q,
                                                        t_QuantityGradients &_grad,
                                                        t_VA_Type *_i_grad)
      {
        __DOASSERT( t_QuantityTraits::size(_q) != coefs.size(),
                       "Wrong number of objective functions during evaluation \n"
                    << "Expected " << t_QuantityTraits::size(_q) << "\n"
                    << "Found " << coefs.size() << "\n" )
        
        t_ScalarQuantity results = 0.0;
        typename t_Quantity :: const_iterator i_val = _q.begin();
        typename t_Quantity :: const_iterator i_val_end = _q.end();
        typename t_QuantityGradients :: iterator i_grad = _grad.begin();
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.end();

        t_VA_Type *i2 = _i_grad + _q.size();
        t_VA_Type *i1 = _i_grad;
        for( ; i1 != i2; ++i1 ) *i1 = t_VA_Type(0);

        t_VA_Type *const_grad_result = new t_VA_Type[ _q.size() ];
        t_VA_Type *const_grad_result_end = const_grad_result + _q.size(); 

        for(; i_objective != i_end and i_val != i_val_end;
              ++i_objective, ++i_coef, ++i_val, ++i_grad )
        {
          i1 = const_grad_result;

          for( ; i1 != const_grad_result_end; ++i1 ) *i1 = t_VA_Type(0);

          results +=   (*i_objective)->evaluate_with_gradient( *i_val, *i_grad, const_grad_result )
                     * (*i_coef);
          i1 = const_grad_result;

          for(; i1 != const_grad_result_end; ++i1, ++i2, ++i_coef ) *i2 +=  (*i_coef) * (*i1);
        }

        delete[] const_grad_result;
          
        return results;
      };
  template< class T_GA_TRAITS >
    void
      LinearSum<T_GA_TRAITS> ::  evaluate_gradient( const t_Quantity &_q,
                                                                t_QuantityGradients &_grad,
                                                                t_VA_Type *_i_grad)
      {
        __DOASSERT( t_QuantityTraits::size(_q) != coefs.size(),
                       "Wrong number of objective functions during evaluation \n"
                    << "Expected " << t_QuantityTraits::size(_q) << "\n"
                    << "Found " << coefs.size() << "\n" )
        
        typename t_Quantity :: const_iterator i_val = _q.begin();
        typename t_Quantity :: const_iterator i_val_end = _q.end();
        typename t_QuantityGradients :: iterator i_grad = _grad.begin();
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.begin();

        t_VA_Type *i2 = _i_grad + _q.size();
        t_VA_Type *i1 = _i_grad;
        for( ; i1 != i2; ++i1 ) *i1 = t_VA_Type(0);

        t_VA_Type *const_grad_result = new t_VA_Type[ _q.size() ];
        t_VA_Type *const_grad_result_end = const_grad_result + _q.size(); 

        for(; i_objective != i_end and i_val != i_val_end;
              ++i_objective, ++i_coef, ++i_val, ++i_grad )
        {
          i1 = const_grad_result;

          for( ; i1 != const_grad_result_end; ++i1 ) *i1 = t_VA_Type(0);

          (*i_objective)->evaluate_gradient( *i_val, *i_grad, const_grad_result );
          i1 = const_grad_result;

          for(; i1 != const_grad_result_end; ++i1, ++i2, ++i_coef ) *i2 +=  (*i_coef) * (*i1);
        }

        delete[] const_grad_result;
      };
  template< class T_GA_TRAITS >
    typename LinearSum<T_GA_TRAITS>::t_VA_Type
      LinearSum<T_GA_TRAITS> :: evaluate_one_gradient( const t_Quantity & _q,
                                                                   t_QuantityGradients& _grad,
                                                                   types::t_unsigned _n)
      {
        __DOASSERT( t_QuantityTraits::size(_q) != coefs.size(),
                       "Wrong number of objective functions during evaluation \n"
                    << "Expected " << t_QuantityTraits::size(_q) << "\n"
                    << "Found " << coefs.size() << "\n" )

        typename t_Quantity :: const_iterator i_val = _q.begin();
        typename t_Quantity :: const_iterator i_val_end = _q.end();
        typename t_QuantityGradients :: iterator i_grad = _grad.begin();
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.begin();

        t_VA_Type result = t_VA_Type(0);
        for(; i_objective != i_end and i_val != i_val_end;
              ++i_objective, ++i_coef, ++i_val, ++i_grad )
          result +=   ( *i_coef ) 
                    * (*i_objective)->evaluate_one_gradient( *i_val, *i_grad, _n );

        return result;
      };
     







  template< class T_GA_TRAITS >
     typename Types<T_GA_TRAITS> :: Scalar*
      Types<T_GA_TRAITS> :: scalar_from_xml( const TiXmlElement &_node )
      {
        if ( not &_node ) return NULL;
        std::string str = "minimize"; 
        std::string name = Print::lowercase(_node.Value());
        if (    name.compare("objective") == 0 
             or name.compare("method") == 0 )
        {
          if ( _node.Attribute( "type" ) )
            str = Print::lowercase(_node.Attribute( "type" ));
        }
        else if ( _node.Attribute("objective") )
          str = _node.Attribute( "objective" );
        if ( str.compare("convexhull") == 0 )
        {
          Print::xmg << Print::Xmg::comment << "Objective: ConvexHull" << Print::endl;
          return new ConvexHull<t_GATraits>;
        }
        else if ( str.compare("minimize") == 0 )
        {
          Print::xmg << Print::Xmg::comment << "Objective: Minimize" << Print::endl;
          return new Minimize<t_GATraits>;
        }
        else if ( str.compare("maximize") == 0 )
        {
          Print::xmg << Print::Xmg::comment << "Objective: Maximize" << Print::endl;
          return new Maximize<t_GATraits>;
        }
        else if (str.compare("target") == 0 )
        {
          if( _node.Attribute("target") )
          {
            double d; _node.Attribute("target", &d );
            Print::xmg << Print::Xmg::comment
                       << "Objective: Target (" << d << ")" << Print::endl;
            return new Target<t_GATraits>( (types::t_real) d );
          }
        }
        if ( _node.FirstChildElement( "Objective" ) )
         return scalar_from_xml( *_node.FirstChildElement( "Objective" ) ); 

        return NULL;
      }
    template< class T_GA_TRAITS >
       typename Types<T_GA_TRAITS> :: Vector*
        Types<T_GA_TRAITS> :: vector_from_xml( const TiXmlElement &_node )
        {
          if ( not &_node ) return NULL;
//         std::string str = "minimize"; 
//         std::string name = Print::lowercase(_node.Value());
//         if (    name.compare("objective") == 0 
//              or name.compare("method") == 0 )
//         {
//           if ( _node.Attribute( "type" ) )
//             str = Print::lowercase(_node.Attribute( "type" ));
//         }
          if ( Vector::t_QuantityTraits::is_vector ) // and str.compare("LinearSum") == 0 )
          {
            LinearSum<T_GA_TRAITS> *linear = new LinearSum<T_GA_TRAITS>;
            __DOASSERT( not linear, "Mememory allocation error\n" )
            
            Print::xmg << Print::Xmg::comment << "Objective: begin LinearSum" << Print::endl;
            Print::xmg << Print::Xmg::indent;
            const TiXmlElement *child = _node.FirstChildElement("Objective");
            for(; child; child = child->NextSiblingElement("Objective") )
            {
              Scalar* scalar = scalar_from_xml( *child );
              if ( not scalar ) continue;
              double d = 0.0;
              if ( not child->Attribute("coef", &d) ) d = 1.0;
              linear->add( scalar, t_ScalarQuantity(d) );
              std::ostringstream sstr; sstr << ", coef=" << d;
              Print::xmg.add_to_last(sstr.str());
            }
            Print::xmg << Print::Xmg::unindent;
            Print::xmg << Print::Xmg::comment << "Objective: end LinearSum" << Print::endl;
            return linear;
          }

          if ( _node.FirstChildElement( "Objective" ) )
           return new_from_xml( *_node.FirstChildElement( "Objective" ) ); 

          return NULL;
        }

} // namespace Objective

} // namspace GA
#endif
