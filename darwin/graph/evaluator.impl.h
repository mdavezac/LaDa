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
          typename Farmer<T_GATRAITS>::t_Individual*
            Farmer<T_GATRAITS> :: onWait( types::t_unsigned _bull )
            {
              typename t_Unknowns :: iterator i_unknown = unknowns.begin();
              typename t_Unknowns :: iterator i_unknown_end = unknowns.end();
              typename t_ProcessIds :: iterator i_id = process_ids.end();
              for(; i_unknown != i_unknown_end; ++i_unknown, ++i_id )
                if( *i_id == _bull )
                {
                  unknowns.erase( i_unknown );
                  process_ids.erase( i_id );
                  break;
                }
              i_unknown = unknowns.begin();
              i_id = process_ids.begin();
              for(; i_unknown != i_unknown_end; ++i_unknown, ++i_id )
                if( *i_id == -1 ) break;
           
              if( i_unknown == i_unknown_end ) return NULL;
              *i_id = _bull;
              return *i_unknown;
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
