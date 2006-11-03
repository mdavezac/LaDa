#ifndef _COLONIZE_H_
#define _COLONIZE_H_

#include <eo/eoEvalFunc.h>
#include <eo/eoFunctor.h>
#include <eo/utils/eoRNG.h>

#include <opt/types.h>
using namespace types;

#include "operators.h"
#include "breeder.h"

namespace LaDa 
{
  template<class t_Object, class t_Islands = std::list< eoPop<t_Object> > > 
  class Colonize : public eoUF<t_Islands, void >
  {
    protected:
      t_unsigned every; // goes to work every "every" generations
      bool stable_pop;  // stable population vs grow population
      eoEvalFunc<t_Object> &eval;
      Breeder<t_Object> &breeder;
      UtterRandom <t_Object> utter_random;

    public:
      Colonize   ( eoEvalFunc<t_Object> &_eval, eoBreed<t_Object> &_breed,
                   t_unsigned _every, bool _stable = false ) 
               : eval(_eval), breeder(_breed), every(_every), stable_pop(_stable) {};
      Colonize   ( const Colonize<t_Object, t_Islands> &_copy )
               : eval(_copy.eval), breeder(_copy.breeder), 
                 every(_copy.every), stable_pop(_copy.stable_pop) {};

      virtual ~Colonize ();

      void operator()( t_Islands &_islands )
      {
        t_unsigned pop_size = _islands.begin()->size();
        t_unsigned new_size = pop_size;
        t_unsigned nb_islands = _islands.size();
        t_unsigned take_from_each;
        t_real d;
        ++nb_islands; 

        if ( stable_pop )
        {
          d = floor( (double) pop_size * (1.0 -  1.0 / ( (double) nb_islands ) ) );
          new_size = ( d > 0 ) ? static_cast<t_unsigned>(d) : 1;
          if ( new_size * nb_islands < pop_size * (nb_islands - 1 ) )
            ++new_size; // won't go to zero
        }
        d = floor( (double) pop_size / (double) ( nb_islands )  ) / 2.0;
        take_from_each = ( d > 0) ? static_cast<t_unsigned>(d) : 1;
        if ( take_from_each * nb_islands < 1  ) 
          take_from_each = 1;
        if ( take_from_each * nb_islands > new_size  )
        {
          new_size = take_from_each * nb_islands;
          d = static_cast<t_real>(new_size) * 0.25 > 2;
          if ( d  < 2 )
            d = 2;
          new_size += static_cast<t_unsigned>(d);
        }

        // swtiches eoHowMany object in breeder.
        eoHowMany howmany(take_from_each);
        eoHowMany **howmany_address = breeder.get_howmany_address();
        eoHowMany *save_howmany = *howmany_address;
        *howmany_address = &howmany;

        // breeds new pop from all other islands
        eoPop<t_Object> new_pop; new_pop.reserve( new_size ); 
        eoPop<t_Object> offsprings; offsprings.reserve( take_from_each );
        typename t_Islands :: iterator i_pop = _islands.begin();
        typename t_Islands :: iterator i_end = _islands.end();
        new_pop.resize( new_size >> 1 );
        typename eoPop<t_Object> :: iterator i_indiv = new_pop.begin();
        for( ; i_pop != i_end; ++i_pop, i_indiv += take_from_each )
        {
          breeder(*i_pop, offsprings);
          std::copy( offsprings.begin(), offsprings.end(), i_indiv );
          offsprings.clear();
        }

        // Now evaluates the new island
        // note that new UtterRandom object are evaluated when created
        // in the loop below
        typename eoPop<t_Object> :: iterator i_indiv_end = new_pop.end();
        for( ; i_indiv != i_indiv_end; ++i_indiv)
          eval(*i_indiv);

        // Shoves island into island container
        _islands.push_back(new_pop);

        // completes the job by resizing each island to proper value
        i_pop = _islands.begin();
        i_end = _islands.end();
        for( ; i_pop != i_end; ++i_pop )
        {
          t_unsigned s = i_pop->size();
          if ( s > new_size )
          {
            i_pop->nth_element(new_size);
            i_pop->resize(new_size);
          }
          else if ( s < new_size )
          {
            t_Object object;
            for(; s < new_size; ++s )
            {
              utter_random(object);
              eval(object); // for simplicity, new objects are evaluated here
              i_pop->push_back(object);
            }
          }
        }

        // and finally, reavaluates baseline stuff if necessary
        if ( not t_Object :: is_baseline_valid() )
        {
          i_pop = _islands.begin();
          for( ; i_pop != i_end; ++i_pop )
          {
            i_indiv = new_pop.begin();
            i_indiv_end = new_pop.end();
            for( ; i_indiv != i_indiv_end; ++i_indiv )
              eval(*i_indiv);
          }
        }

        // returns breeder to what it was
        *howmany_address =  save_howmany;

      }

  };

} // namespace LaDa

#endif
