#ifndef _CHECKPOINT_H_
#define _CHECKPOINT_H_

#include <eo/utils/eoUpdater.h>

#include "taboo.h"

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
  template< class t_Object >
  class NuclearWinter : public eoStatBase<t_Object>
  {
    protected:
      Taboo_Base<t_Object> &taboo;
      TriggeredOp<t_Object> &triggered_op;
      unsigned nuclear_winter_length, nuclear_winter_age;
      bool is_gone_nuclear;

    public:
      NuclearWinter   ( Taboo_Base<t_Object> &_taboo, 
                        TriggeredOp<t_Object> &_triggered_op,
                        unsigned _nwl )
                    : taboo(_taboo),
                      triggered_op(_triggered_op),
                      nuclear_winter_length(_nwl),
                      is_gone_nuclear(false) {}

      virtual ~NuclearWinter(){};

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
            triggered_op.trigger(false); 
            is_gone_nuclear = false;
            std::cout << "Nuclear-winter is over" << std::endl;
          }
        }
        else if ( taboo.is_problematic() )
        { // starts nuclear winter
          nuclear_winter_age = 0;
          is_gone_nuclear = true;
          triggered_op.trigger(true); 
          taboo.set_problematic(false);
          std::cout << "Going Nuclear" << std::endl;
        }
      }


      // some anoying stuff
      void printOn( std::ostream &__os ) const {};
      void readFrom( std::istream &__os ) const {};
      void lastCall() {};
      virtual std::string className(void) const { return "LaDa::NuclearWinter"; }

  };

  template<class t_Object>
  class UpdateAgeTaboo : public eoStatBase<t_Object> 
  {
    protected:
      Taboo< t_Object, std::list<t_Object> > & taboo;
      eoIncrementorParam<unsigned> &age;
      unsigned max_age;

    public:
      UpdateAgeTaboo  ( Taboo< t_Object, std::list<t_Object> > & _taboo,
                        eoIncrementorParam<unsigned> &_age,
                        unsigned _max_age )
                     : taboo(_taboo), age(_age), max_age(_max_age) {};
      UpdateAgeTaboo  ( const UpdateAgeTaboo<t_Object> & _update )
                     : taboo(_update.taboo), age(_update.age),
                       max_age(_update.max_age) {};

      virtual void operator()( const eoPop<t_Object> &_pop )
      {
        typename eoPop<t_Object> :: const_iterator i_pop = _pop.begin();
        typename eoPop<t_Object> :: const_iterator i_end = _pop.end();

        for( ; i_pop != i_end; ++i_pop )
          if ( age.value() - i_pop->get_age() > max_age )
            taboo.add( *i_pop );
      }

      // some anoying stuff
      void printOn( std::ostream &__os ) const {};
      void readFrom( std::istream &__os ) const {};
      void lastCall()
        { taboo.print_out( std::cout ); }

      virtual std::string className(void) const { return "LaDa::UpdateTaboo"; }
  };

} // namespace LaDa


#endif //  _CHECKPOINT_H_
