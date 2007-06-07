#ifndef _POPGROWTH_H_
#define _POPGROWTH_H_

#include <eo/eoEvalFunc.h>
#include <eo/eoFunctor.h>
#include <eo/utils/eoRNG.h>

#include <opt/types.h>
using namespace types;
#include "eotypes.h"

#include "operators.h"
#include "breeder.h"

namespace LaDa 
{
  template<class t_Object, class t_Islands = std::list< eoPop<t_Object> > > 
  class PopGrowth : public eoUF<t_Islands&, void >
  {
    protected:
      t_unsigned every; // goes to work every "every" generations
      t_unsigned growth_rate;
      t_unsigned max_pop;
      eoEvalFunc<t_Object> &eval;
      Breeder<t_Object> &breeder;
      t_unsigned nb_count;
      eoPop<t_Object> &offsprings;

    public:
      PopGrowth   ( eoEvalFunc<t_Object> &_eval, Breeder<t_Object> &_breed,
                   t_unsigned _every, t_unsigned _gr, t_unsigned _max_pop,
                   eoPop<t_Object> &_offsprings ) 
                : every(_every), growth_rate(_gr), max_pop(_max_pop), eval(_eval),
                  breeder(_breed), nb_count(0), offsprings(_offsprings) {};
      PopGrowth   ( const Colonize<t_Object, t_Islands> &_copy )
                : every(_copy.every), growth_rate(_copy.growth_rate),
                  max_pop(_copy.max_pop), eval(_copy.eval), breeder(_copy.breeder), 
                  nb_count(_copy.nb_count), offsprings(_copy.offsprings) {};

      virtual ~PopGrowth () {};

      void operator()( t_Islands &_islands )
      {
        ++nb_count;
        if ( nb_count % every )
          return;

        // swtiches eoHowMany object in breeder.
        eoHowMany **howmany_address = breeder.get_howmany_address();
        eoHowMany *save_howmany = *howmany_address;
        eoHowMany howmany(static_cast<eotypes::t_int>(growth_rate));
        *howmany_address = &howmany;

        // breeds new pop from all other islands
        typename t_Islands :: iterator i_pop = _islands.begin();
        typename t_Islands :: iterator i_end = _islands.end();
        for( ; i_pop != i_end; ++i_pop )
          if ( i_pop->size() < max_pop )
          {
            // creates new objects
            offsprings.clear();
            breeder(*i_pop, offsprings);

            // makes sure we don't add too many objects
            t_unsigned osize = offsprings.size();
            t_unsigned psize = i_pop->size();
            if ( osize + psize > max_pop )
              offsprings.resize( max_pop - psize );

            // evaluates and appends new objects
            i_pop->reserve( i_pop->size() + offsprings.size() );
            typename eoPop<t_Object> :: iterator i_indiv = offsprings.begin();
            typename eoPop<t_Object> :: iterator i_indiv_end = offsprings.end();
            for(; i_indiv != i_indiv_end; ++i_indiv)
            {
              eval(*i_indiv);
              i_pop->push_back( *i_indiv );
            }
          }

        // and finally, reavaluates baseline stuff if necessary
        if ( not t_Object :: is_baseline_valid() )
        {
          i_pop = _islands.begin();
          for( ; i_pop != i_end; ++i_pop )
          {
            typename eoPop<t_Object> :: iterator i_indiv = i_pop->begin();
            typename eoPop<t_Object> :: iterator i_indiv_end = i_pop->end();
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
