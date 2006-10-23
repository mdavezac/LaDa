#ifndef _CHECKPOINT_H_
#define _CHECKPOINT_H_

#include <eo/eoOp.h>
#include <eo/eoGenOp.h>
#include <eo/utils/eoUpdater.h>
#include <eo/utils/eoHowMany.h>
#include <eo/utils/eoStat.h>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "taboo.h"

#include<opt/types.h>
using namespace types;

namespace LaDa 
{
  // implements a call back which simply calls a print_xmgrace and a print_xml
  template<class t_Call_Back>
  class PrintXmgrace : public eoUpdater
  {
    private:
      t_Call_Back* call_back;
    public:
      PrintXmgrace( t_Call_Back* _call_back) { call_back =_call_back; }
      virtual void lastCall()
        { call_back->print_xml(); }
      virtual void operator()(void)
        { call_back->print_xmgrace(); }
      
      virtual std::string className(void) const { return "Monitor"; } 
  };

  // checks for taboo unconvergence from a breeder, 
  // response. If response does not get through 
  template< class t_Object, class t_Call_Back >
  class NuclearWinter : public eoStatBase<t_Object>
  {
    protected:
      Taboo_Base<t_Object> &taboo;
      eoGenOp<t_Object> &normal_ops;
      eoSequentialOp<t_Object> nuclear_ops;
      eoGenOp<t_Object> **breeding_ops;
      t_unsigned nuclear_winter_length, nuclear_winter_age;
      t_Call_Back &call_back;
      bool is_gone_nuclear;
      eoHowMany nuclear_howmany;
      eoHowMany normal_howmany;
      eoHowMany **breeding_howmany;

    public:
      NuclearWinter   ( Taboo_Base<t_Object> &_taboo, 
                        eoGenOp<t_Object> &_nops,
                        eoGenOp<t_Object> &_nwops,
                        t_Call_Back &_call_back,
                        t_real &_normal_howmany )
                    : taboo(_taboo),
                      normal_ops(_nops),
                      breeding_ops(NULL),
                      nuclear_winter_length(2),
                      call_back( _call_back ),
                      is_gone_nuclear(false),
                      nuclear_howmany(1.0),
                      normal_howmany(_normal_howmany),
                      breeding_howmany(NULL)
      {
        nuclear_ops.add(_nops, 1.0);
        nuclear_ops.add(_nwops, 1.0);
      }

      virtual ~NuclearWinter(){};

      void set_howmany( eoHowMany **_howmany)
      {
        breeding_howmany = _howmany;
        *breeding_howmany = &normal_howmany;
      }

      // checks for taboo unconvergenced
      // if true, then transforms each individual according to response
      // if still cannot create new individuals, then returns false
      // and GA exits
      virtual void operator()( const eoPop<t_Object> &_pop )
      {
        if ( is_gone_nuclear )
        { // ends nuclear winter
          ++nuclear_winter_age;
          if ( nuclear_winter_age > nuclear_winter_length ) 
          {
            *breeding_ops = &normal_ops; 
            *breeding_howmany = &normal_howmany;
            is_gone_nuclear = false;
            std::string str("Nuclear-winter is over");
            call_back.print_xmgrace(str);
          }
        }
        else if ( taboo.is_problematic() )
        { // starts nuclear winter
          nuclear_winter_age = 0;
          is_gone_nuclear = true;
          *breeding_ops = &nuclear_ops; 
          *breeding_howmany = &nuclear_howmany;
          taboo.set_problematic(false);
          std::string str("Going Nuclear");
          call_back.print_xmgrace(str);
        }
      }

      void set_op_address( eoGenOp<t_Object> ** _ops )
        { breeding_ops = _ops; } 

      eoGenOp<t_Object>* get_op_address() const
      {
        if (is_gone_nuclear)
          return &normal_ops;
        return &nuclear_ops;
      }

      // some anoying stuff
      void printOn( std::ostream &__os ) const {};
      void readFrom( std::istream &__os ) const {};
      virtual void lastCall( const eoPop<t_Object> &_pop) {};
      virtual std::string className(void) const { return "LaDa::NuclearWinter"; }

  };

  template<class t_Object, class t_Call_Back>
  class UpdateAgeTaboo : public eoStatBase<t_Object> 
  {
    protected:
      Taboo< t_Object, std::list<t_Object> > & taboo;
      eoIncrementorParam<t_unsigned> &age;
      t_Call_Back &call_back;
      t_unsigned max_age;
      t_unsigned check_every;
      bool do_print_out;

    public:
      UpdateAgeTaboo  ( Taboo< t_Object, std::list<t_Object> > & _taboo,
                        eoIncrementorParam<t_unsigned> &_age, t_Call_Back &_call_back,
                        t_unsigned _max_age, bool _do_print_out = false )
                     : taboo(_taboo), age(_age), call_back(_call_back), max_age(_max_age),
                       do_print_out( _do_print_out )
      {
        check_every = max_age / 10;
        if ( check_every == 0 )
          check_every = 1;
      };
      UpdateAgeTaboo  ( const UpdateAgeTaboo<t_Object, t_Call_Back> & _update )
                     : taboo(_update.taboo), age(_update.age), call_back( _update.call_back ),
                       max_age(_update.max_age), check_every( _update.check_every ),
                       do_print_out(_update.do_print_out) {};

      virtual void operator()( const eoPop<t_Object> &_pop )
      {
        t_unsigned ga_age = age.value();
        if ( ga_age < max_age or ga_age%check_every != 0 )
          return;

        typename eoPop<t_Object> :: const_iterator i_pop = _pop.begin();
        typename eoPop<t_Object> :: const_iterator i_end = _pop.end();

        for( ; i_pop != i_end; ++i_pop )
          if ( ga_age - i_pop->get_age() > max_age )
          {
            taboo.add( *i_pop );
            if ( do_print_out )
            { 
              std::ostringstream sstr;
              sstr << " AgeTaboo: new individual ";
              i_pop->print_out( sstr );
              std::string str = sstr.str();
              call_back.print_xmgrace( str );
            }
          }
      }

      // some anoying stuff
      void printOn( std::ostream &__os ) const {};
      void readFrom( std::istream &__os ) const {};
      virtual void lastCall( const eoPop<t_Object> &_pop)
        {} //;taboo.print_out( std::cout ); }

      virtual std::string className(void) const { return "LaDa::UpdateTaboo"; }
  };

  // Gets an average of individuals accumulated over all generations
  template< class t_Object, class t_Call_Back >
  class AccAverage : public eoStatBase<t_Object>
  {
    protected:
      t_Call_Back &call_back;
      eoIncrementorParam<t_unsigned> &age;
      t_unsigned total_individuals;
      t_Object average;

    public:
      AccAverage   ( t_Call_Back &_call_back,
                     eoIncrementorParam<t_unsigned> &_age,
                     t_unsigned _size )
                 : call_back( _call_back ), age(_age),
                   total_individuals(0)
      {           
        average.resize( _size ); // should initialize to 0 in opt_function_base.h
      }
      virtual ~AccAverage() {}
      virtual std::string className(void) const { return "LaDa::AccAverage"; }
      virtual void operator()( const eoPop<t_Object> &_pop )
      {
        t_unsigned this_age = age.value();
        typename eoPop<t_Object> :: const_iterator i_indiv = _pop.begin();
        typename eoPop<t_Object> :: const_iterator i_end = _pop.end();
        typename t_Object :: iterator i_av_begin = average.begin();
        typename t_Object :: iterator i_av_end = average.end();
        typename t_Object :: iterator i_average;
        typename t_Object :: const_iterator i_var;
        for( ; i_indiv != i_end; ++i_indiv )
          if ( i_indiv->get_age() == this_age )
          {
            ++total_individuals;
            i_var = i_indiv->begin();
            i_average = i_av_begin ;
            for ( ; i_average != i_av_end; ++i_var, ++i_average )
               *i_average += *i_var;
          }

        // prints stuff out
        {
          std::ostringstream sstr; 
          typename t_Object :: const_iterator i_var = average.begin();
          typename t_Object :: const_iterator i_var_end = average.end();
          sstr << "AccAverage: " << std::setw(5) << std::setprecision(2);
          for( ; i_var != i_var_end; ++i_var )
            sstr << (*i_var / (t_real) total_individuals) << " ";
          std::string str = sstr.str();
          call_back.print_xmgrace( str );
        }
      }
  };

  // Gets an average of individuals accumulated over this generations
  template< class t_Object, class t_Call_Back >
  class PopAverage : public eoStatBase<t_Object>
  {
    protected:
      t_Call_Back &call_back;
      t_Object average;

    public:
      PopAverage   ( t_Call_Back &_call_back,
                     t_unsigned _size )
                 : call_back( _call_back )
      {
        average.resize( _size ); // should initialize to 0 in opt_function_base.h
      }
      virtual ~PopAverage() {}
      virtual std::string className(void) const { return "LaDa::PopAverage"; }
      virtual void operator()( const eoPop<t_Object> &_pop )
      {
        t_unsigned pSize = _pop.size();
        typename t_Object :: iterator i_av_begin = average.begin();
        typename t_Object :: iterator i_av_end = average.end();
        typename t_Object :: iterator i_average;
        typename eoPop<t_Object> :: const_iterator i_indiv = _pop.begin();
        typename eoPop<t_Object> :: const_iterator i_end = _pop.end();
        typename t_Object :: const_iterator i_var;
        average = *i_indiv;
        for( ++i_indiv; i_indiv != i_end; ++i_indiv )
        {
          i_average = i_av_begin;
          i_var = i_indiv->begin();
          for ( ; i_average != i_av_end; ++i_var, ++i_average )
             *i_average += *i_var;
        }

        // prints stuff out
        {
          std::ostringstream sstr; 
          typename t_Object :: const_iterator i_var = average.begin();
          typename t_Object :: const_iterator i_var_end = average.end();
          sstr << "PopAverage: " << std::setw(5) << std::setprecision(2);
          for( ; i_var != i_var_end; ++i_var )
            sstr << (*i_var / (t_real) pSize ) << " ";
          std::string str = sstr.str();
          call_back.print_xmgrace( str );
        }
      }
  };

  // island continuator. Takes two population iterators and goes through them
  template<class t_Object, class t_Islands = std::list< eoPop<t_Object> > > 
  class IslandsContinuator : public eoContinue<t_Object>
  {
    typedef typename t_Islands :: iterator iterator;
    typedef typename t_Islands :: const_iterator const_iterator;
    protected:
      std::list < eoContinue<t_Object>* >       continuators;
      std::list < eoSortedStatBase<t_Object>* > sorted;
      std::list < eoStatBase<t_Object>* >       stats;
      std::list < eoMonitor* >                  monitors;
      std::list < eoUpdater* >                  updaters;

    public:
      IslandsContinuator( eoContinue<t_Object> &_cont ) 
        { continuators.push_back( &_cont ); }

      void add(eoContinue<t_Object>& _cont)
        { continuators.push_back(&_cont); }
      void add(eoSortedStatBase<t_Object>& _stat)
        { sorted.push_back(&_stat); }
      void add(eoStatBase<t_Object>& _stat) 
        { stats.push_back(&_stat); }
      void add(eoMonitor& _mon)        
        { monitors.push_back(&_mon); }
      void add(eoUpdater& _upd)        
        { updaters.push_back(&_upd); }

      bool operator()(const eoPop<t_Object>& _pop)
      { return true; }

      // applies population dependant stuff
      void apply_stats(const eoPop<t_Object>& _pop) const
      {
        // first goes through sorted stats
        if ( not sorted.empty() )
        { 
          typename std::list < eoSortedStatBase<t_Object>* > :: const_iterator i_sorted = sorted.begin();
          typename std::list < eoSortedStatBase<t_Object>* > :: const_iterator i_end = sorted.end();
          std::vector<const t_Object*> sorted_pop;
          _pop.sort(sorted_pop);
          for ( ; i_sorted != i_end; ++i_sorted )
            (*i_sorted)->operator()(sorted_pop);
        }
        // then through unsorted stats
        if ( not stats.empty() )
        { 
          typename std::list < eoStatBase<t_Object>* > :: const_iterator i_stat = stats.begin();
          typename std::list < eoStatBase<t_Object>* > :: const_iterator i_end = stats.end();
          for ( ; i_stat != i_end; ++i_stat )
            (*i_stat)->operator()(_pop);
        }
      } // bool apply ( ... )

      void apply_monitors_updaters()
      {
        // applies monitors
        if ( not monitors.empty() )
        { 
          typename std::list < eoMonitor* > :: const_iterator i_monitor = monitors.begin();
          typename std::list < eoMonitor* > :: const_iterator i_end = monitors.end();
          for ( ; i_monitor != i_end; ++i_monitor )
            (*i_monitor)->operator()();
        }
        // then through updaters
        if ( not updaters.empty() )
        { 
          typename std::list < eoUpdater* > :: const_iterator i_updater = updaters.begin();
          typename std::list < eoUpdater* > :: const_iterator i_end = updaters.end();
          for ( ; i_updater != i_end; ++i_updater )
            (*i_updater)->operator()();
        }
      }

      bool apply_continuators( const eoPop<t_Object> & _pop ) const
      {
        typename std::list < eoContinue<t_Object>* > :: const_iterator i_continue = continuators.begin();
        typename std::list < eoContinue<t_Object>* > :: const_iterator i_end = continuators.end();
        bool result = (*i_continue)->operator()( _pop); // generation counter
        for ( ++i_continue; i_continue != i_end; ++i_continue )
          if ( not (*i_continue)->operator()( _pop) )
            result = false;
        return result;
      }

      void lastCall( const eoPop<t_Object>& _pop ) const 
      {
        // first goes through sorted stats
        if ( not sorted.empty() )
        { 
          typename std::list < eoSortedStatBase<t_Object>* > :: const_iterator i_sorted = sorted.begin();
          typename std::list < eoSortedStatBase<t_Object>* > :: const_iterator i_end = sorted.end();
          std::vector<const t_Object*> sorted_pop;
          _pop.sort(sorted_pop);
          for ( ; i_sorted != i_end; ++i_sorted )
            (*i_sorted)->lastCall(sorted_pop);
        }
        // then through unsorted stats
        if ( not stats.empty() )
        { 
          typename std::list < eoStatBase<t_Object>* > :: const_iterator i_stat = stats.begin();
          typename std::list < eoStatBase<t_Object>* > :: const_iterator i_end = stats.end();
          for ( ; i_stat != i_end; ++i_stat )
            (*i_stat)->lastCall(_pop);
        }
        // applies monitors
        if ( not monitors.empty() )
        { 
          typename std::list < eoMonitor* > :: const_iterator i_monitor = monitors.begin();
          typename std::list < eoMonitor* > :: const_iterator i_end = monitors.end();
          for ( ; i_monitor != i_end; ++i_monitor )
            (*i_monitor)->lastCall();
        }
        // then through updaters
        if ( not updaters.empty() )
        { 
          typename std::list < eoUpdater* > :: const_iterator i_updater = updaters.begin();
          typename std::list < eoUpdater* > :: const_iterator i_end = updaters.end();
          for ( ; i_updater != i_end; ++i_updater )
            (*i_updater)->lastCall();
        }
      }

      bool apply ( iterator &_i_begin, iterator &_i_end )
      {
        iterator i_pop = _i_begin;
        bool result = true;

        for( ; i_pop != _i_end; ++i_pop )
          apply_stats( *i_pop );

        apply_monitors_updaters();

        for( i_pop = _i_begin; i_pop != _i_end; ++i_pop )
          if ( not apply_continuators( *i_pop ) )
            result = false;
 
        if ( not result )
          for( i_pop = _i_begin; i_pop != _i_end; ++i_pop )
            lastCall( *i_pop );

        return result;
      }

      virtual std::string className(void) const { return "LaDa::IslandsContinuator"; }
  };

} // namespace LaDa


#endif //  _CHECKPOINT_H_
