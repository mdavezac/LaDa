//
//  Version: $Id$
//
#ifndef _DARWIN_STATISTICS_H_
#define _DARWIN_STATISTICS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <eo/eoPop.h>
#include <eo/eoGenOp.h>

#include "opt/types.h"
#include "print/xmg.h"

namespace GA
{
  template< class T_INDIVIDUAL>
  class TrueCensus : public eoStatBase< T_INDIVIDUAL >
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
    protected:
      typedef eoPop<t_Individual> t_Population;
    private:
      typedef typename t_Individual :: t_Object t_Object;
      
    public:
      TrueCensus () {}
      virtual ~TrueCensus() {}
      virtual std::string className(void) const { return "GA::TrueCensus"; }
      virtual void operator()( const t_Population &_pop )
      {
        typename t_Population :: const_iterator i_begin = _pop.begin();
        typename t_Population :: const_iterator i_indiv1 = i_begin;
        typename t_Population :: const_iterator i_end = _pop.end();
        typename t_Population :: const_iterator i_indiv2;
        types::t_unsigned N = _pop.size();
        
        for( ; i_indiv1 != i_end; ++i_indiv1 )
        {
          i_indiv2 = i_indiv1; ++i_indiv2;
          if ( i_end != std::find_if( i_indiv2, i_end,
                                      std::bind1st(std::equal_to<t_Object>(), *i_indiv1) ) )
            --N;
        }

        // prints stuff out
        Print::xmg << Print::Xmg::comment <<  "True Census: "
                   << N << " / " << _pop.size() << Print::endl; 
      }
  };
  
}

#endif
