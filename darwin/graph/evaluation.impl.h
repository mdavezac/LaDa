//
//  Version: $Id$
//

#include <boost/tuple/tuple.hpp>
#include <boost/serialization/utility.hpp>

namespace GA
{
  namespace mpi 
  {
    namespace Graph
    {
      namespace Evaluation
      {
        template<class T_GATRAITS, template< class > class T_BASE >
          void Farmer<T_GATRAITS, T_BASE> :: operator()( t_Population &_pop,
                                                         t_Population &_offspring ) 
            { 
              __ASSERT(     t_Base :: objective
                        and t_CommBase :: objective, 
                        "objective pointer has not been set.\n" )
              // eventually, this will call Evaluator::Farmer::evaluate(), 
              // and put unknown individuals into a list.
              t_Base::evaluate(_offspring); 


              // At this point,  unknowns are evaluated. 
              t_CommBase :: start_all();
              while( notdone() ) t_CommBase::test_bulls();
              t_CommBase :: wait_bulls(); 
              t_CommBase::comm->barrier();

              // Ok, now recomputes everything through base class. 
              // if invalid, recomputes whole population
              t_Base::operator()( _pop, _offspring );
            }
        template<class T_GATRAITS, template< class > class T_BASE >
          typename Farmer<T_GATRAITS, T_BASE> :: t_FitnessQuantity
            Farmer<T_GATRAITS, T_BASE> :: evaluate( t_Individual &_indiv )
            {
              if ( not _indiv.invalid() ) return t_Base::evaluate( _indiv );
              
              // fitness AND quantities of _indiv must be valid from here-on
              ++t_Base::nb_eval;
              __DODEBUGCODE( Print::out << "Evaluating Individual: " << _indiv << Print::endl; )
              unknowns.push_back( t_Unknown( &_indiv, -1) );
              return t_FitnessQuantity(0);
            }

        template<class T_GATRAITS, template< class > class T_BASE >
          void Farmer<T_GATRAITS, T_BASE> :: onWait( types::t_unsigned _bull )
          {
            t_Individual *indiv = next( _bull );
            if( not indiv )
            {
              t_CommBase::send_command( _bull, t_CommBase::t_Commands::DONE );
              return;
            }
            t_CommBase::send_command( _bull, t_CommBase::t_Commands::GO );
            t_CommBase::comm->send( _bull, ONWAIT_TAG( TAG ), *indiv );
            t_CommBase::activate( _bull );
          }

        template<class T_GATRAITS, template< class > class T_BASE >
          typename Farmer<T_GATRAITS, T_BASE> :: t_Individual*
            Farmer<T_GATRAITS, T_BASE> :: next( types::t_unsigned _bull )
            {
              namespace blamb = boost::lambda;
              typename t_Unknowns :: iterator i_unknown
                 = std::find_if( unknowns.begin(), unknowns.end(), 
                                   blamb::bind( &t_Unknown::second, blamb::_1 )
                                == blamb::constant( _bull ) );
              if ( i_unknown != unknowns.end() )
              {
                typedef ::mpi::Pair
                        < 
                          typename t_GATraits :: t_QuantityTraits :: t_Quantity&,
                          typename t_GATraits :: t_IndivTraits :: t_Object &
                        > t_MPIPair;
                t_MPIPair mpipair( i_unknown->first->quantities(),
                                   i_unknown->first->Object() );
                t_CommBase::comm->recv( _bull, ONWAIT_TAG( TAG+1 ), mpipair );
//               t_CommBase::comm->recv
//               (
//                 _bull, ONWAIT_TAG( TAG+1 ), 
//                 i_unknown->first->quantities()
//               );
//               t_CommBase::comm->recv
//               (
//                 _bull, ONWAIT_TAG( TAG+2 ), 
//                 i_unknown->first->Object() 
//               );
                __TRYDEBUGCODE(
                    t_Base::objective->init( *i_unknown->first );
                    i_unknown->first->set_fitness(
                     (*t_Base::objective)( i_unknown->first->const_quantities() ) );,
                    "Error while evaluating fitness.\n" )
                unknowns.erase( i_unknown );
              }
  
              i_unknown = std::find_if( unknowns.begin(), unknowns.end(), 
                                           blamb::bind( &t_Unknown::second, blamb::_1 )
                                        == blamb::constant(-1) );
              if( i_unknown == unknowns.end() ) return NULL;
              i_unknown->second = _bull;
              return i_unknown->first;
            }



        template<class T_GATRAITS, template< class > class T_BASE >
          void Bull<T_GATRAITS, T_BASE> :: operator()(t_Population& _parents,
                                                      t_Population& _offspring)
          {
            t_CommBase :: request( t_CommBase::t_Requests::WAITING );
            while ( t_CommBase::obey() != t_CommBase :: t_Commands :: DONE )
            {
              t_Individual individual;
              t_CommBase::comm->recv( 0, ONWAIT_TAG( TAG ), individual );
              metaeval.init( individual );
              metaeval.evaluate();
              typedef ::mpi::Pair
                      < 
                        typename t_GATraits :: t_QuantityTraits :: t_Quantity&,
                        typename t_GATraits :: t_IndivTraits :: t_Object &
                      > t_MPIPair;
              t_MPIPair mpipair( individual.quantities(),
                                 individual.Object() );
              t_CommBase::comm->send( 0, ONWAIT_TAG( TAG+1 ), mpipair );
//             t_CommBase::comm->send
//             ( 
//               0, ONWAIT_TAG( TAG+1 ),
//               individual.quantities()
//             ); 
//             t_CommBase::comm->send
//             ( 
//               0, ONWAIT_TAG( TAG+2 ),
//               individual.Object() 
//             ); 
              t_CommBase :: request( t_CommBase::t_Requests::WAITING );
            }
            t_CommBase::command( t_CommBase::t_CowCommands::DONE );
            t_CommBase::comm->barrier();
          }

      } // namespace Evaluation
    } // namespace Graph
  } // namespace mpi
} // namespace GA
