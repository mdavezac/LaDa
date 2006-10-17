//-----------------------------------------------------------------------------

#ifndef _BREEDER_H_
#define _BREEDER_H_

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

/**
  Base class for breeders using generalized operators.
*/

namespace LaDa 
{
  template<class t_Object>
  class Breeder: public eoBreed<t_Object>
  {
    protected:
      eoSelectOne<t_Object>& select;
      eoGenOp<t_Object> *op; 
      eoIncrementorParam<unsigned> &age;
      eoHowMany howMany;
    
    public:
      Breeder   ( eoSelectOne<t_Object>& _select, eoGenOp<t_Object>& _op, 
                  eoIncrementorParam<unsigned> &_age,
                  double  _rate=1.0, bool _interpret_as_rate = true) 
              : select( _select ), op(&_op), age(_age),
                howMany(_rate, _interpret_as_rate) {}
      Breeder   ( eoSelectOne<t_Object>& _select, eoGenOp<t_Object>& _op,
                  eoIncrementorParam<unsigned> &_age, eoHowMany &_howmany )
              : select( _select ), op(&_op),
                age(_age), howMany(_howmany) {}
      Breeder   ( Breeder<t_Object> & _breeder )
              : select( _breeder.select ), op(_breeder.op),
                age(_breeder.age), howMany(_breeder.howmany) {}

      virtual ~Breeder() {};
     
      void operator()(const eoPop<t_Object>& _parents, eoPop<t_Object>& _offspring)
      {
        unsigned target = howMany(_parents.size());
     
        _offspring.clear();
        eoSelectivePopulator<t_Object> it(_parents, _offspring, select);
     
        while (_offspring.size() < target)
        {
          (*op)(it);
          (*it).set_age(age.value());
          ++it;
        }
     
        _offspring.resize(target);   // you might have generated a few more
      }

      eoGenOp<t_Object>** get_op_address() 
        { return &op; }
     
      /// The class name.
      virtual std::string className() const { return "LaDa::Breeder"; }

  };

} // namespace Breeder

#endif

