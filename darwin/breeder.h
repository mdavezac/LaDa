//-----------------------------------------------------------------------------

#ifndef _DARWIN_BREEDER_H_
#define _DARWIN_BREEDER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//-----------------------------------------------------------------------------

/*****************************************************************************
 * eoGeneralBreeder: transforms a population using the generalOp construct.
 *****************************************************************************/

#include <eo/eoOp.h>
#include <eo/eoGenOp.h>
#include <eo/eoPopulator.h>
#include <eo/eoSelectOne.h>
#include <eo/eoBreed.h>
#include <eo/utils/eoHowMany.h>

#include "opt/types.h"

#include "checkpoints.h"
/**
  Base class for breeders using generalized operators.
*/

namespace darwin 
{
  template<class T_INDIVIDUAL, class T_POPULATION = eoPop<T_INDIVIDUAL> >
  class Breeder: public eoBreed<T_INDIVIDUAL>
  {
    protected:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_POPULATION t_Population;

    protected:
      eoSelectOne<t_Individual>& select;
      eoGenOp<t_Individual> *op; 
      GenCount &age;
      eoHowMany *howMany;
      eoHowMany *howMany_save;
    
    public:
      Breeder   ( eoSelectOne<t_Individual>& _select, eoGenOp<t_Individual>& _op,
                  GenCount &_age )
              : select( _select ), op(&_op),
                age(_age), howMany(NULL), howMany_save(NULL) {}
      Breeder   ( Breeder<t_Individual> & _breeder )
              : select( _breeder.select ), op(_breeder.op),
                age(_breeder.age), howMany(_breeder.howMany), howMany_save(NULL) {}

      virtual ~Breeder() { if( howMany_save ) delete howMany_save; };
     
      void operator()(const t_Population& _parents, t_Population& _offspring)
      {
        types::t_unsigned target = (*howMany)( (types::t_unsigned)_parents.size());
#ifdef _MPI
        types::t_int residual = target % (mpi::main.size());
        target = target / mpi::main.size();
        if ( mpi::main.rank() < residual )
          ++target;
#endif
     
        _offspring.clear();
        eoSelectivePopulator<t_Individual> it(_parents, _offspring, select);
     
        while (_offspring.size() < target)
        {
          (*op)(it);
          (*it).set_age(age());
          ++it;
        }
     
        _offspring.resize(target);   // you might have generated a few more
      }
#ifdef _MPI // "All Gather" offsprings -- note that replacement schema is deterministic!!
      void synchronize_offsprings( t_Population &_pop )
      {
        mpi::AllGather allgather( mpi::main );
        typename t_Population::iterator i_indiv = _pop.begin();
        typename t_Population::iterator i_indiv_end = _pop.end();
        for (; i_indiv != i_indiv_end; ++i_indiv )
          i_indiv->broadcast(allgather);
        allgather.allocate_buffers();
        for (i_indiv = _pop.begin(); i_indiv != i_indiv_end; ++i_indiv )
          i_indiv->broadcast(allgather);
        allgather();
        _pop.clear();
        t_Individual indiv;
        while( indiv.broadcast( allgather ) )
          { _pop.push_back(indiv); }
      }
#endif 

      eoGenOp<t_Individual>** get_op_address() 
        { return &op; }
      eoHowMany** get_howmany_address() 
        { return &howMany; }
      void set_howmany(types::t_real _rate)
      {
        if ( howMany_save ) delete howMany_save;
        howMany_save = new eoHowMany(_rate);
        howMany = howMany_save;
      }
     
      /// The class name.
      virtual std::string className() const { return "Darwin::Breeder"; }
  };

} // namespace Breeder

#endif

