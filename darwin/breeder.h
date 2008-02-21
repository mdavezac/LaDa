//
//  Version: $Id$
//
//-----------------------------------------------------------------------------

#ifndef _DARWIN_BREEDER_H_
#define _DARWIN_BREEDER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//-----------------------------------------------------------------------------

#include <eo/eoOp.h>
#include <eo/eoGenOp.h>
#include <eo/eoPopulator.h>
#include <eo/eoSelectOne.h>
#include <eo/eoBreed.h>
#include <eo/utils/eoHowMany.h>

#include <opt/types.h>
#include <mpi/mpi_object.h>

#include "checkpoints.h"
#include "gatraits.h"
#include "debug.h"

namespace GA 
{
  //! \brief Creates an offpsring population from a current population
  //! \details Breeder simply accretes all data necessary to the creation of
  //! offspring via mating. It expects to be given a mating operator on input,
  //! which it will then call until the requested number of offspring is
  //! created. It also expects an eoSelectOne object on input for selecting
  //! parents.
  template<class T_GATRAITS>
  class Breeder: public eoBreed<typename T_GATRAITS::t_Individual>
  {
    public:
      typedef T_GATRAITS t_GATraits; //!< all %GA traits
    protected:
      //!  all individual traits
      typedef typename t_GATraits::t_IndivTraits t_IndivTraits;
       //!  type of an individual
      typedef typename t_GATraits::t_Individual  t_Individual; 
        //! type of the population
      typedef typename t_GATraits::t_Population  t_Population;

    protected:
      //! A selection operator for the obtention of parents
      eoSelectOne<t_Individual>* select;
      eoGenOp<t_Individual> *op;  //!< Mating operator
      GenCount *age;              //!< Generation counter
      //! \brief Number of offspring to create
      //! \details It is declared a pointer so that it can be changed
      //! dynamically during the run. This is useful mostly for
      //! GA::NuclearWinter.
      eoHowMany *howMany;     
      //! \brief place holder 
      eoHowMany *howMany_save;
    
    public:
      //! Constructor and Initializer
      Breeder   ( eoSelectOne<t_Individual>* _select,
                  eoGenOp<t_Individual>* _op,
                  GenCount *_age )
              : select( _select ), op(_op),
                age(_age), howMany(NULL), howMany_save(NULL) {}
      //! Copy Constructor
      Breeder   ( Breeder<t_Individual>& _breeder )
              : select( _breeder.select ), op(_breeder.op),
                age(_breeder.age), howMany(_breeder.howMany), howMany_save(NULL) {}
      //! Destructor
      virtual ~Breeder() { if( howMany_save ) delete howMany_save; };
     
      //! Creates \a _offspring population from \a _parent
      void operator()(const t_Population& _parents, t_Population& _offspring);


      //! \brief Hack, don't use
      //! \details Used to change the mating operators from right under the nose of
      //! a Breeeder instance. Used by GA::NuclearWinter.
      //! \todo remove this dirty hack
      eoGenOp<t_Individual>** get_op_address()  { return &op; }
      //! \brief Hack, don't use
      //! \details Used to change the number of offspring to create from right
      //! under the nose of  a Breeeder instance. Used by GA::NuclearWinter.
      //! \todo remove this dirty hack
      eoHowMany** get_howmany_address()  { return &howMany; }
      //! Sets the selector.
      void set( eoSelectOne<t_Individual> *_select ){ select = _select; }
      //! Sets the breeding operators
      void set( eoSelectOne<t_Individual> *_op ){ op = _op; }
      //! Sets the replacement rate
      void set( types::t_real _rep );
      //! Sets the generation counter.
      void set( GenCount *_age ){ age = _age; }
     
      ///! The class name. EO required
      virtual std::string className() const { return "Darwin::Breeder"; }
  };

  template<class T_GATRAITS>
  inline void Breeder<T_GATRAITS> :: operator()(const t_Population& _parents,
                                                t_Population& _offspring)
  {
    types::t_unsigned target = (*howMany)( (types::t_unsigned) _parents.size());
    _offspring.clear();
    eoSelectivePopulator<t_Individual> it(_parents, _offspring, *select);
  
    while (_offspring.size() < target)
    {
      (*op)(it);
      (*it).set_age( (*age)() );
      ++it;
    }
  
    _offspring.resize(target);   // you might have generated a few more
  }
  template<class T_GATRAITS>
  inline void Breeder<T_GATRAITS> :: set(types::t_real _rate)
  {
    if ( howMany_save ) delete howMany_save;
    howMany_save = new eoHowMany(_rate);
    howMany = howMany_save;
  }

} // namespace Breeder

#endif

