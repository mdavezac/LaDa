//
//  Version: $Id$
//

#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>


namespace LaDa
{
  namespace GA
  {
    namespace mpi 
    {
      namespace Graph
      {
        namespace Evaluator
        {
          template<class T_GATRAITS>
            inline void Bull<T_GATRAITS> :: evaluate()
            {
              __ASSERT( not trueEvaluator, "True Evaluator is not set." )
              t_CommBase :: command( t_CommBase::t_CowCommands::EVALUATE );
              boost::mpi::broadcast( *t_CommBase :: cowcomm,
                                     *t_Base::current_individual, 0 );

              trueEvaluator->evaluate();
            }
          template<class T_GATRAITS>
            inline void Bull<T_GATRAITS>
              :: evaluate_gradient( t_QuantityGradients& _grad )
            {
              __ASSERT( not trueEvaluator, "True Evaluator is not set." )
              t_CommBase :: command( t_CommBase::t_CowCommands::GRADIENT );
              boost::mpi::broadcast( *t_CommBase :: cowcomm,
                                     *t_Base::current_individual, 0 );
              trueEvaluator->evaluate_gradient( _grad );
            }
          template<class T_GATRAITS>
            inline void Bull<T_GATRAITS>
              :: evaluate_with_gradient( t_QuantityGradients& _grad )
            {
              __ASSERT( not trueEvaluator, "True Evaluator is not set." )
              t_CommBase :: command( t_CommBase::t_CowCommands::WITH_GRADIENT );
              boost::mpi::broadcast( *t_CommBase :: cowcomm,
                                     *t_Base::current_individual, 0 );
              trueEvaluator->evaluate_with_gradient( _grad );
            }
          template<class T_GATRAITS>
            inline void Bull<T_GATRAITS>
              :: evaluate_one_gradient( t_QuantityGradients& _grad,
                                        types::t_unsigned _pos )
            {
              __ASSERT( not trueEvaluator, "True Evaluator is not set." )
              t_CommBase :: command( t_CommBase::t_CowCommands::ONE_GRADIENT );
              boost::mpi::broadcast( *t_CommBase :: cowcomm,
                                     *t_Base::current_individual, 0 );
              boost::mpi::broadcast( *t_CommBase :: cowcomm, _pos, 0 );
              trueEvaluator->evaluate_one_gradient( _grad, _pos );
            }

        } // namespace Evaluation
      } // namespace Graph
    } // namespace mpi
  } // namespace GA
} // namespace LaDa
