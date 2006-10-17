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
      eoGenOp<t_Object> &normal_ops;
      eoSequentialOp<t_Object> nuclear_ops;
      eoGenOp<t_Object> **breeding_ops;
      unsigned nuclear_winter_length, nuclear_winter_age;
      bool is_gone_nuclear;

    public:
      NuclearWinter   ( Taboo_Base<t_Object> &_taboo, 
                        eoGenOp<t_Object> &_nops,
                        eoGenOp<t_Object> &_nwops,
                        unsigned _nwl )
                    : taboo(_taboo),
                      normal_ops(_nops),
                      breeding_ops(NULL),
                      nuclear_winter_length(_nwl),
                      is_gone_nuclear(false)
      {
        nuclear_ops.add(_nops, 1.0);
        nuclear_ops.add(_nwops, 1.0);
      }

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
            *breeding_ops = &normal_ops; 
            is_gone_nuclear = false;
            std::cout << "Nuclear-winter is over" << std::endl;
          }
        }
        else if ( taboo.is_problematic() )
        { // starts nuclear winter
          nuclear_winter_age = 0;
          is_gone_nuclear = true;
          *breeding_ops = &nuclear_ops; 
          taboo.set_problematic(false);
          std::cout << "Going Nuclear" << std::endl;
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

  template<class t_Object>
  class UpdateAgeTaboo : public eoStatBase<t_Object> 
  {
    protected:
      Taboo< t_Object, std::list<t_Object> > & taboo;
      eoIncrementorParam<unsigned> &age;
      unsigned max_age;
      unsigned check_every;

    public:
      UpdateAgeTaboo  ( Taboo< t_Object, std::list<t_Object> > & _taboo,
                        eoIncrementorParam<unsigned> &_age,
                        unsigned _max_age )
                     : taboo(_taboo), age(_age), max_age(_max_age) 
      {
        check_every = max_age / 10;
        if ( check_every == 0 )
          check_every = 1;
      };
      UpdateAgeTaboo  ( const UpdateAgeTaboo<t_Object> & _update )
                     : taboo(_update.taboo), age(_update.age),
                       max_age(_update.max_age), check_every( _update.check_every ) {};

      virtual void operator()( const eoPop<t_Object> &_pop )
      {
        unsigned ga_age = age.value();
        if ( ga_age < max_age or ga_age%check_every != 0 )
          return;

        typename eoPop<t_Object> :: const_iterator i_pop = _pop.begin();
        typename eoPop<t_Object> :: const_iterator i_end = _pop.end();

        for( ; i_pop != i_end; ++i_pop )
          if ( ga_age - i_pop->get_age() > max_age )
            taboo.add( *i_pop );
      }

      // some anoying stuff
      void printOn( std::ostream &__os ) const {};
      void readFrom( std::istream &__os ) const {};
      virtual void lastCall( const eoPop<t_Object> &_pop)
        {  } //;taboo.print_out( std::cout ); }

      virtual std::string className(void) const { return "LaDa::UpdateTaboo"; }
  };

} // namespace LaDa


#endif //  _CHECKPOINT_H_
