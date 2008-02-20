//
//  Version: $Id$
//


#include <list>

namespace GA
{
  namespace mpi 
  {
    namespace Graph
    {
      //! \brief Contains all evalaution related stuff in the mpi::Graph Topology.
      //! \details The classes below are meant to be used as derived from
      //!          "standard" Evaluation classes. They will simply overwrite
      //!          the behavior of the Evaluation::Base::operator()(
      //!          t_Population &_parent, t_Population &_offspring ) routine.
      //!          Instead of the standard approach to evaluating offsprings,
      //!          each invalid individual is dispatched by the farmer to the
      //!          herds for evaluation. Individuals are dispatched as herds
      //!          become available. The classes from which the following
      //!          classes are derive is defined through a template. These
      //!          classes also derive from the Comm classes, \e via the CRT
      //!          mechanism.
      namespace Evaluation
      {
        template<class T_BASE>
        void Farmer<T_BASE> :: operator()(const t_Population& _parents,
                                          t_Population& _offspring)
        {
          __ASSERT( evaluator, "Pointer to Evaluator is not set.\n")
          __ASSERT( objective, "Pointer to Evaluator is not set.\n")
          __ASSERT( store, "Pointer to Evaluator is not set.\n")

          std::list<t_Individual &> unknowns;
          typename t_Population :: iterator i_indiv = _offspring.begin();
          typename t_Population :: iterator i_indiv_end = _offspring.end();
          for(; i_indiv != i_indiv_end; ++i_indiv )
            if( i_indiv->invalid() ) unknowns.push_back( i_indiv );

          while( unknowns.size() ) t_CommBase::test_bulls();
        }

      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA
