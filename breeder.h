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

#include <opt/types.h>
using namespace types;
#include <eo/eotypes.h>

#include "checkpoint.h"
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
      GenCount &age;
      eoHowMany *howMany;
      eoHowMany *howMany_save;
    
    public:
      Breeder   ( eoSelectOne<t_Object>& _select, eoGenOp<t_Object>& _op,
                  GenCount &_age )
              : select( _select ), op(&_op),
                age(_age), howMany(NULL), howMany_save(NULL) {}
      Breeder   ( Breeder<t_Object> & _breeder )
              : select( _breeder.select ), op(_breeder.op),
                age(_breeder.age), howMany(_breeder.howMany), howMany_save(NULL) {}

      virtual ~Breeder() { if( howMany_save ) delete howMany_save; };
     
      void operator()(const eoPop<t_Object>& _parents, eoPop<t_Object>& _offspring)
      {
        t_unsigned target = (*howMany)( (eotypes::t_unsigned)_parents.size());
     
        _offspring.clear();
        eoSelectivePopulator<t_Object> it(_parents, _offspring, select);
     
        while (_offspring.size() < target)
        {
          (*op)(it);
          (*it).set_age(age());
          ++it;
        }
     
        _offspring.resize(target);   // you might have generated a few more
      }

      eoGenOp<t_Object>** get_op_address() 
        { return &op; }
      eoHowMany** get_howmany_address() 
        { return &howMany; }
      void set_howmany(t_real _rate)
      {
        if ( howMany_save ) delete howMany_save;
        howMany_save = new eoHowMany(_rate);
        howMany = howMany_save;
      }
     
      /// The class name.
      virtual std::string className() const { return "LaDa::Breeder"; }

  };

} // namespace Breeder

#endif

