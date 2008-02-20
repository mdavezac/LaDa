//
//  Version: $Id$
//



namespace GA
{
  namespace mpi 
  {
    namespace Graph
    {
      namespace Evaluator
      {
        template<class T_GATRAITS>
          typename Farmer<T_GATRAITS>::t_Individual
            Farmer<T_GATRAITS> :: onWait( types::t_unsigned _bull )
            {
              typename t_Unknowns :: iterator i_unknown = unknowns.begin();
              typename t_Unknowns :: iterator i_unknown_end = unknowns.end();
              for(; i_unknown != i_unknown_end; ++i_unknown )
                if( i_unknown->second == _bull - 1 )
                {
                  unkowns.erase( i_unknown );
                  break;
                }
              i_unknown = unknowns.begin();
              for(; i_unknown != i_unknown_end; ++i_unknown )
                if( i_unknown->second == -1 ) break;
           
              if( i_unknown == i_unknown_end ) return NULL;
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
          inline void Bull<T_GATRAITS> :: evaluate_gradient( t_QuantityGradients& _grad )
          {
            t_CommBase :: command( t_CommBase::t_CowCommands::GRADIENT );
            t_CommBase :: bcast( *t_Base::current_individual );
            t_Base::evaluate_gradient();
          }
        template<class T_GATRAITS>
          inline void Bull<T_GATRAITS> :: evaluate_with_gradient( t_QuantityGradients& _grad )
          {
            t_CommBase :: command( t_CommBase::t_CowCommands::WITH_GRADIENT );
            t_CommBase :: bcast( *t_Base::current_individual );
            t_Base::evaluate_with_gradient( _grad, _pos );
          }
        template<class T_GATRAITS>
          inline void Bull<T_GATRAITS> :: evaluate_one_gradient( t_QuantityGradients& _grad,
                                                                 types::t_unsigned &_pos )
          {
            t_CommBase :: command( t_CommBase::t_CowCommands::ONE_GRADIENT );
            t_CommBase :: bcast( *t_Base::current_individual );
            t_CommBase :: bcast( pos );
            t_Base::evaluate_one_gradient( _grad, _pos );
          }

      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA
