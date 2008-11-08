//
//  Version: $Id$
//

namespace LaDa
{
  namespace GA
  {
    namespace mpi 
    {
      namespace Graph
      {
          template< class T_GA_TRAITS >
            const typename BullObjective<T_GA_TRAITS>::t_Fitness&
              BullObjective<T_GA_TRAITS> :: operator()( const t_Quantity& _q ) 
              {
                t_CommBase :: request( t_CommBase :: t_Requests :: OBJECTIVE );
                t_CommBase::comm->send( 0, ONOBJECTIVE_TAG1( TAG ), 
                                        *t_Base::current_indiv );
                t_CommBase::comm->recv( 0, ONOBJECTIVE_TAG2( TAG ), t_Base::fitness );
                return t_Base::fitness;
              };

          template< class T_GA_TRAITS >
            typename BullObjective<T_GA_TRAITS>::t_ScalarQuantity
              BullObjective<T_GA_TRAITS>
                :: evaluate_with_gradient( const t_Quantity &_q,
                                           t_QuantityGradients &_grad,
                                           t_VA_Type *_i_grad)
              {
                t_CommBase::comm->send( 0, ONWITHGRADIENT_TAG1( TAG ), 
                                        *t_Base::current_indiv );
                t_CommBase::comm->send( 0, ONWITHGRADIENT_TAG2( TAG ), _grad );
                t_CommBase::comm->recv( 0, ONWITHGRADIENT_TAG3( TAG ),
                                        _i_grad, _grad.size() );
                types::t_real result;
                t_CommBase::comm->recv( 0, ONWITHGRADIENT_TAG4( TAG ), result );
                return result;
              };

          template< class T_GA_TRAITS >
            void BullObjective<T_GA_TRAITS>
              :: evaluate_gradient( const t_Quantity &_q,
                                    t_QuantityGradients &_grad,
                                    t_VA_Type *_i_grad)
            {
              t_CommBase::comm->send( 0, ONGRADIENT_TAG1( TAG ), 
                                      *t_Base::current_indiv );
              t_CommBase::comm->send( 0, ONGRADIENT_TAG2( TAG ), _grad );
              t_CommBase::comm->recv( 0, ONGRADIENT_TAG3( TAG ), _i_grad, _grad.size() );
            };

          template< class T_GA_TRAITS >
            typename BullObjective<T_GA_TRAITS>::t_VA_Type
              BullObjective<T_GA_TRAITS>
                :: evaluate_one_gradient( const t_Quantity &_q,
                                          t_QuantityGradients &_grad,
                                          types::t_unsigned _n)
              {
                t_CommBase::comm->send( 0, ONONEGRADIENT_TAG1( TAG ), 
                                        *t_Base::current_indiv );
                t_CommBase::comm->send( 0, ONONEGRADIENT_TAG2( TAG ), _grad );
                t_CommBase::comm->send( 0, ONONEGRADIENT_TAG3( TAG ), _n );
                t_VA_Type result;
                t_CommBase::comm->recv( 0, ONONEGRADIENT_TAG4( TAG ), result );
                return result;
              };
      } // namespace Graph
    } // namespace mpi
  } // namespace GA
} // namespace LaDa
