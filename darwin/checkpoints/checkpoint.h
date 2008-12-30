//
//  Version: $Id: checkpoints.h 849 2008-11-12 04:45:34Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_CHECKPOINT_H_
#define _LADA_GA_CHECKPOINTS_CHECKPOINT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <boost/signal.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/function_traits.hpp>
#include <boost/static_assert.hpp>
#include <boost/function_types/function_arity.hpp>
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>



#include <mpi/mpi_object.h>
#include <opt/debug.h>

void shit() {};

namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      namespace details
      {
        //! Returns  "... and valueN and valueN+1 ..." on signal return values.
        struct BoolSignalsNor
        {
          //! Return type.
          typedef bool result_type;

          template<typename T_ITERATOR>
          bool operator()(T_ITERATOR _first, T_ITERATOR _last) const
          {
            // nothing in signal, throw errror.
            __ASSERT( _first == _last, "No continuator. "
                     "The GA doesn't know when to stop.\n" )
      
            bool result = *_first++;
            for(; _first != _last; ++_first )
              result &= *_first;
          
            return result;
          }
        };
      }

      //! \brief Aggregates and calls all checkpoint types.
      //! \details The three checkpoint type differ by their arguments and functions:
      //!            . continuators are of type bool(). If any one continuator
      //!              returns false, the GA should stop.
      //!            . updaters are of type void(bool) where the arguments if
      //!              this is the last call (a continuator returned false).
      //!            . statistics are of type void(bool, const t_Population&).
      //!              First argument specifies whether this is the last call.
      template< class T_ISLANDS >
        class CheckPoint
        {
          public:
            //! Type of the islands.
            typedef T_ISLANDS t_Islands;
            //! Type of the population.
            typedef typename t_Islands :: value_type t_Population;
            //! Type of the continuator functors.
            typedef bool (t_Continuator)(void);
            //! Type of the updater functors.
            typedef void (t_Updater)(bool);
            //! Type of the statistic functors.
            typedef void (t_Statistic)(bool, const t_Population&);

          public:

            //! Constructor.
            CheckPoint() {}
            //! Destructor.
            virtual ~CheckPoint() {}

            //! Calls continuators, updaters, statistics.
            bool operator()( const t_Islands& _islands );

            //! Connects updater functors.
            template< class T_FUNCTOR > 
              void connect_updater( const T_FUNCTOR& _functor )
                { updaters_.connect( _functor ); }
            //! Connects continuators functors.
            template< class T_FUNCTOR > 
              void connect_continuator( const T_FUNCTOR& _functor )
                { continuators_.connect( _functor ); }
            //! Connects continuators functors.
            template< class T_FUNCTOR > 
              void connect_statistic( const T_FUNCTOR& _functor )
                { statistics_.connect( _functor ); }

            //! Connects age counter. 
            void connect_age_counter( const GenCount& _age ) { age_ = _age; }

          protected:
            //! Type of the continuator container.
            typedef boost::signal< t_Continuator, details::BoolSignalsNor > t_Continuators;
            //! Type of the updater container.
            typedef boost::signal< t_Updater > t_Updaters;
            //! Type of the update container.
            typedef boost::signal< t_Statistic > t_Statistics;

            //! The continuators container.
            t_Continuators continuators_;
            //! The updaters container.
            t_Updaters updaters_;
            //! The statistics container.
            t_Statistics statistics_;
            //! Reference to generational counter.
            GenCount age_;
        };


      template< class T_ISLANDS >
        bool CheckPoint<T_ISLANDS> :: operator()( const t_Islands& _islands )
        {
          bool keepgoing( continuators_() );
          __MPICODE
          (
            boost::mpi::communicator world;
            keepgoing = boost::mpi::all_reduce( world, keepgoing, std::logical_and<bool>() ); 
          )
          const bool lastcall( not keepgoing );
          updaters_( lastcall );
          typename t_Islands :: const_iterator i_pop = _islands.begin();
          typename t_Islands :: const_iterator i_pop_end = _islands.end();
          for(; i_pop != i_pop_end; ++i_pop )
            statistics_( lastcall, *i_pop );
          ++age_;
          return keepgoing;
        } 

    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif
