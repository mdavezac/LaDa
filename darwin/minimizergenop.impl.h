//
//  Version: $Id$
//
#ifndef _MINIMIZERGENOP_IMPL_H_
#define _MINIMIZERGENOP_IMPL_H_

namespace LaDa
{
  namespace GA
  {

    template<class T_GATRAITS>
      inline typename Minimizer_Functional<T_GATRAITS> :: t_VA_Type 
        Minimizer_Functional<T_GATRAITS> :: evaluate()
        { 
          evaluation->evaluate( *current_indiv ); 
          return (t_VA_Type) current_indiv->fitness(); 
        }
    template<class T_GATRAITS>
      inline typename Minimizer_Functional<T_GATRAITS> :: t_VA_Type 
        Minimizer_Functional<T_GATRAITS> :: evaluate_with_gradient( t_VA_Type *_i_grad )
        {
          evaluation->evaluate_with_gradient( *current_indiv,
                                              gradients, _i_grad ); 
          return current_indiv->fitness(); 
        }
    template<class T_GATRAITS>
      inline typename Minimizer_Functional<T_GATRAITS> :: t_VA_Type 
        Minimizer_Functional<T_GATRAITS> :: evaluate_one_gradient( types::t_unsigned _pos )
        {
          types::t_real result = evaluation->evaluate_one_gradient( *current_indiv, gradients, _pos); 
          return result;
        }
    template<class T_GATRAITS>
      inline bool Minimizer_Functional<T_GATRAITS> :: init( t_Individual & _indiv)
      {
        savestate.init(_indiv);
        evaluation->init(_indiv); 
        variables = &_indiv.Object().Container();
        gradients.resize( variables->size(), _indiv.const_quantities() );
        Traits::zero_out( gradients );
        current_indiv = &_indiv;
        return true;
      }


    template< class T_GATRAITS >
      inline void SaveStateIndividual<T_GATRAITS> :: save()
      {
        if ( not current_indiv ) return;
        quantity = current_indiv->const_quantities();
        fitness =  current_indiv->fitness();
      }
    template< class T_GATRAITS >
      inline void SaveStateIndividual<T_GATRAITS> :: reset() const 
      {
        if ( not current_indiv ) return;
        current_indiv->quantities() = quantity;
        current_indiv->fitness() = fitness;
      }



    template< class T_GATRAITS >
      inline void MinimizerGenOp<T_GATRAITS> :: apply(eoPopulator<t_Individual>& _pop)
      {
        functional.init( *_pop );
        (*minimizer)( functional );
        functional.evaluate();
      }
    template< class T_GATRAITS>
    bool MinimizerGenOp<T_GATRAITS> :: Load( const TiXmlElement &_node )
    {
      std::string name = _node.Value();
      if ( name.compare("Minimizer") != 0 )
        return false;
     
      if ( not _node.Attribute("type") )
        return false;
      
      name = _node.Attribute("type");
     
      if ( name.compare("VA") == 0 )
      {
        Print::xmg << Print::Xmg::comment << "VA optimizer" << Print::endl;
        // pointer is owned by caller !!
        // don't deallocate
        minimizer =  new Minimizer::VA<t_Functional, t_SaveState>( _node, functional.savestate );
      }
      else if ( name.compare("SA") == 0 )
      {
        Print::xmg << Print::Xmg::comment << "SA optimizer" << Print::endl;
        // pointer is owned by caller !!
        // don't deallocate
        minimizer = new Minimizer::VA<t_Functional, t_SaveState>( _node, functional.savestate );
      }
      else if ( name.compare("Beratan") == 0 )
      {
        Print::xmg << Print::Xmg::comment << "Beratan optimizer" << Print::endl;
        // pointer is owned by caller !!
        // don't deallocate
        minimizer = new Minimizer::Beratan<t_Functional>( _node );
      }
     
      return minimizer != NULL;
     
    }


  } // namespace GA
} // namespace LaDa

#endif // _MINIMIZERGENOP_IMPL_H_
