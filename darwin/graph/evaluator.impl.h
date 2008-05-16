//
//  Version: $Id$
//

#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>


namespace GA
{
  namespace mpi 
  {
    namespace Graph
    {
      namespace Evaluator
      {
        template<class T_GATRAITS>
          typename Farmer<T_GATRAITS>::t_Individual*
            Farmer<T_GATRAITS> :: onWait( types::t_unsigned _bull )
            {
              namespace blamb = boost::lambda;
              typename t_Unknowns :: iterator i_unknown
                 = std::find_if( unknowns.begin(), unknowns.end(), 
                                   blamb::bind( &t_Unknown::second, blamb::_1 )
                                == blamb::constant( _bull ) );
              if( i_unknown == unknowns.end() )
                __THROW_ERROR( "Individual evaluation returned from unassigned process" );
              unknowns.erase( i_unknown );


              i_unknown = std::find_if( unknowns.begin(), unknowns.end(), 
                                           blamb::bind( &t_Unknown::second, blamb::_1 )
                                        == blamb::constant(-1) );
              if( i_unknown == unknowns.end() ) return NULL;
              i_unknown->second = _bull;
              return i_unknown->first;
            }
          

        template<class T_GATRAITS>
          inline void Bull<T_GATRAITS> :: evaluate()
          {
            t_CommBase :: command( t_CommBase::t_CowCommands::EVALUATE );
            t_CommBase :: bcast( *t_Base::current_individual );
            t_Base::evaluate();
          }
        template<class T_GATRAITS>
          inline void Bull<T_GATRAITS>
            :: evaluate_gradient( t_QuantityGradients& _grad )
          {
            t_CommBase :: command( t_CommBase::t_CowCommands::GRADIENT );
            t_CommBase :: bcast( *t_Base::current_individual );
            t_Base::evaluate_gradient( _grad );
          }
        template<class T_GATRAITS>
          inline void Bull<T_GATRAITS>
            :: evaluate_with_gradient( t_QuantityGradients& _grad )
          {
            t_CommBase :: command( t_CommBase::t_CowCommands::WITH_GRADIENT );
            t_CommBase :: bcast( *t_Base::current_individual );
            t_Base::evaluate_with_gradient( _grad );
          }
        template<class T_GATRAITS>
          inline void Bull<T_GATRAITS>
            :: evaluate_one_gradient( t_QuantityGradients& _grad,
                                      types::t_unsigned _pos )
          {
            t_CommBase :: command( t_CommBase::t_CowCommands::ONE_GRADIENT );
            t_CommBase :: bcast( *t_Base::current_individual );
            t_CommBase :: bcast( _pos );
            t_Base::evaluate_one_gradient( _grad, _pos );
          }

      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA
