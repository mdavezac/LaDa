//
//  Version: $Id$
//
#ifndef _DARWIN_OPERATORS_H_
#define _DARWIN_OPERATORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <exception>
#include <iostream>
#include <vector>
#include <list>
#include <utility>

#include <eo/eoPop.h>
#include <eo/eoOp.h>
#include <eo/eoOpContainer.h>
#include <eo/utils/eoRNG.h>



#include "opt/types.h"
#include "gencount.h"
#include "debug.h"

namespace GA 
{
  //! \brief Sequential container of genetic operators.
  //! \details Stores pointers to genetic operators. With each operator is
  //!          associated a rate. For each operator, a coined biased according
  //!          to the rate is flipped. The operator is called only if the toss
  //!          returns true.
  //! \pre rates should be in the range [0,1]
  template<class T_GATRAITS>
  class SequentialOp : public eoOpContainer<typename T_GATRAITS :: t_Individual >
  {
    protected:
      //! All relevant %GA types.
      typedef T_GATRAITS t_GATraits;
    protected:
      //! The type of the individual
      typedef typename t_GATraits :: t_Individual t_Individual;
      //! The type of the base class
      typedef eoOpContainer< t_Individual > t_Base;
      using t_Base :: ops;
      using t_Base :: rates;

    public:
      //! Constructor
      SequentialOp() : t_Base() {}

      //! The genetic operator itself
      virtual void apply(eoPopulator<t_Individual> &_populator);
    
      //! Returns "GA::SequentialOp"
      virtual std::string className() const {return "GA::SequentialOp";}

  };

  template<class T_GATRAITS>
    inline void SequentialOp<T_GATRAITS> :: apply(eoPopulator<t_Individual> &_populator)
    {
      Print :: out << "SequentialOp 0 " << Print::endl;
      typename std::vector< types::t_real > :: const_iterator i_rate = rates.begin();
      typename std::vector< eoGenOp<t_Individual>* > :: iterator i_op = ops.begin();
      typename std::vector< eoGenOp<t_Individual>* > :: iterator i_end = ops.end();
      Print :: out << "SequentialOp 1 " << Print::endl;
      for (; i_op != i_end; ++i_op, ++i_rate )
      {
      Print :: out << "SequentialOp 2 " << Print::endl;
        if ( eo::rng.flip( *i_rate ) )
        {
      Print :: out << "SequentialOp 3 " << typeid( *(*i_op) ).name() << Print::endl;
          (*i_op)->operator()( _populator );
      Print :: out << "SequentialOp 4 " << Print::endl;
        }
      Print :: out << "SequentialOp 5 " << Print::endl;
      }
      Print :: out << "SequentialOp 6 " << Print::endl;
    }



  //! \brief Proportional container of genetic operators.
  //! \details Stores pointers to genetic operators. With each operator is
  //!          associated a rate. A roulette is performed with these rates to
  //!          determine which of the genetic operators to called. 
  template<class T_GATRAITS>
  class ProportionalOp : public eoOpContainer<typename T_GATRAITS :: t_Individual >
  {
    protected:
      //! All relevant %GA types.
      typedef T_GATRAITS t_GATraits;
    protected:
      //! The type of the individual
      typedef typename t_GATraits :: t_Individual t_Individual;
      //! The type of the base class
      typedef eoOpContainer< t_Individual > t_Base;
      using t_Base :: ops;
      using t_Base :: rates;

    public:
      //! Constructor
      ProportionalOp() : t_Base() {}
      
      //! The genetic operator itself
      virtual void apply(eoPopulator<t_Individual> &_populator);

      //! Returns "GA::ProportionalOp"
      virtual std::string className() const {return "GA::ProportionalOp";}
  };

  template<class T_GATRAITS>
    inline void ProportionalOp<T_GATRAITS> :: apply(eoPopulator<t_Individual> &_populator)
    {
      types::t_unsigned i = static_cast<types::t_unsigned>(eo::rng.roulette_wheel(rates));
      return (*ops[i])( _populator );
    }



  //! \brief Wraps an on/off switch around an eoOp operator.
  //! \details The operator is only called if TriggeredOp::is_triggered is
  //!          true. This variable can be set by a call to
  //!          TriggeredOp::trigger.
  template<class T_GATRAITS>
  class TriggeredOp : public eoGenOp<typename T_GATRAITS :: t_Individual >
  {
    public:
      //! All relevant %GA types.
      typedef T_GATRAITS t_GATraits;
    protected:
      //! The type of the individual
      typedef typename t_GATraits :: t_Individual t_Individual;

    protected:
      //! True if TriggeredOp::op is to be called.
      bool is_triggered;
      //! Points to  the operator to call when trigger is on.
      eoGenOp<t_Individual> *op;

    public:
      //! \brief Constructor.
      //! \param _op is wrapped into a genetic operator, TriggeredOp::op.
      //! \param _store is used to store a genetic operator wrapped around \a _op
      //! \param _t true if the operator should on at start. False by default.
      TriggeredOp  ( eoOp<t_Individual> &_op,
                     eoFunctorStore &_store,
                     bool _t = false  )
                  : is_triggered(_t)
        { op = &wrap_op<t_Individual>( _op, _store ); }
      //! Destructor
      virtual ~TriggeredOp() {};
    
      //! The maximum number of produced offspring.
      virtual types::t_unsigned max_production()
        { return op->max_production(); }

      //! Switches the operator on or off. (off=default).
      void trigger( bool _trigger = false )
       { is_triggered = _trigger; }
      //! Returns "GA::TriggeredOps"
      virtual std::string className() const {return "GA::TriggeredOps";}

      //! The genetic operator functor.
      virtual void apply( eoPopulator<t_Individual> &_indiv ) 
        { if ( is_triggered ) (*op)(_indiv); }

  };

  //! \brief Calls an eoOp operator periodically.
  //! \details The operator is called once every PeriodicOp::period. More
  //!          specifically, this class references a GA::GenCount object. When
  //!          this object returns a multiple of PeriodicOp::period, then the
  //!          operator to which Periodic::op points is called.
  template<class T_GATRAITS>
  class PeriodicOp : public eoGenOp<typename T_GATRAITS :: t_Individual >
  {
    public:
      //! All relevant %GA types.
      typedef T_GATRAITS t_GATraits;
    protected:
      //! The type of the individuals
      typedef typename t_GATraits :: t_Individual t_Individual;
      
    protected:
      //! Period of the calls
      types::t_unsigned period;
      //! The "calendar" object. 
      GenCount &age;
      //! Points to the operator to call.
      eoGenOp<t_Individual> *op;
    public:
      //! \brief Constructor.
      //! \param _op is wrapped into a genetic operator, PeriodicOp::op.
      //! \param _period is the period with which to call \a _op.
      //! \param _age should reference the generation counter of the run.
      //! \param _store is used to store a genetic operator wrapped around \a _op.
      PeriodicOp   ( eoOp<t_Individual> &_op,
                     types::t_unsigned _period, 
                     GenCount &_age,
                     eoFunctorStore &_store)
                 : period(_period), age(_age)
        { op = &wrap_op<t_Individual>( _op, _store ); }
      //! Destructor
      virtual ~PeriodicOp() {};
    
      //! Returns the number of produced offpsring
      virtual types::t_unsigned max_production()
        { return op->max_production(); }

      //! Returns "GA::PeriodicOp".
      virtual std::string className() const {return "GA::PeriodicOp";}

      //! This is the genetic functor.
      virtual void apply( eoPopulator<t_Individual> &_indiv ) 
        { if ( age() % period == 0 ) (*op)(_indiv); }
  };

} // endif GA

#endif
