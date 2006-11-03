#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#include <eo/utils/eoStat.h>
#include <eo/eoPop.h>

#include <algorithm>
#include <complex>
#include <cmath>

#include<opt/types.h>
using namespace types;

namespace LaDa
{
  // Gets an average of individuals accumulated over all generations
  template< class t_Call_Back, class t_Object = typename t_Call_Back :: t_Object >
  class AccAverage : public eoStatBase<t_Object>
  {
    protected:
      t_Call_Back &call_back;
      GenCount &age;
      t_unsigned total_individuals;
      t_Object average;

    public:
      AccAverage   ( t_Call_Back &_call_back,
                     GenCount &_age,
                     t_unsigned _size )
                 : call_back( _call_back ), age(_age),
                   total_individuals(0)
      {           
        average.resize( _size ); // should initialize to 0 in opt_function_base.h
      }
      virtual ~AccAverage() {}
      virtual std::string className(void) const { return "LaDa::AccAverage"; }
      virtual void operator()( const eoPop<t_Object> &_pop )
      {
        t_unsigned this_age = age();
        typename eoPop<t_Object> :: const_iterator i_indiv = _pop.begin();
        typename eoPop<t_Object> :: const_iterator i_end = _pop.end();
        typename t_Object :: iterator i_av_begin = average.begin();
        typename t_Object :: iterator i_av_end = average.end();
        typename t_Object :: iterator i_average;
        typename t_Object :: const_iterator i_var;
        for( ; i_indiv != i_end; ++i_indiv )
          if ( i_indiv->get_age() == this_age )
          {
            ++total_individuals;
            i_var = i_indiv->begin();
            i_average = i_av_begin ;
            for ( ; i_average != i_av_end; ++i_var, ++i_average )
               *i_average += *i_var;
          }

        // prints stuff out
        {
          std::ostringstream sstr; 
          typename t_Object :: const_iterator i_var = average.begin();
          typename t_Object :: const_iterator i_var_end = average.end();
          sstr << "AccAverage: " << std::setw(5) << std::setprecision(2);
          for( ; i_var != i_var_end; ++i_var )
            sstr << (*i_var / (t_real) total_individuals) << " ";
          std::string str = sstr.str();
          call_back.print_xmgrace( str );
        }
      }
  };

  // Gets an average of individuals accumulated over this generations
  template< class t_Call_Back, class t_Object = typename t_Call_Back :: t_Object >
  class PopAverage : public eoStatBase<t_Object>
  {
    protected:
      t_Call_Back &call_back;
      t_Object average;

    public:
      PopAverage   ( t_Call_Back &_call_back,
                     t_unsigned _size )
                 : call_back( _call_back )
      {
        average.resize( _size ); // should initialize to 0 in opt_function_base.h
      }
      virtual ~PopAverage() {}
      virtual std::string className(void) const { return "LaDa::PopAverage"; }
      virtual void operator()( const eoPop<t_Object> &_pop )
      {
        t_unsigned pSize = _pop.size();
        typename t_Object :: iterator i_av_begin = average.begin();
        typename t_Object :: iterator i_av_end = average.end();
        typename t_Object :: iterator i_average;
        typename eoPop<t_Object> :: const_iterator i_indiv = _pop.begin();
        typename eoPop<t_Object> :: const_iterator i_end = _pop.end();
        typename t_Object :: const_iterator i_var;
        average = *i_indiv;
        for( ++i_indiv; i_indiv != i_end; ++i_indiv )
        {
          i_average = i_av_begin;
          i_var = i_indiv->begin();
          for ( ; i_average != i_av_end; ++i_var, ++i_average )
             *i_average += *i_var;
        }

        // prints stuff out
        {
          std::ostringstream sstr; 
          typename t_Object :: const_iterator i_var = average.begin();
          typename t_Object :: const_iterator i_var_end = average.end();
          sstr << "PopAverage: " << std::setw(5) << std::setprecision(2);
          for( ; i_var != i_var_end; ++i_var )
            sstr << (*i_var / (t_real) pSize ) << " ";
          std::string str = sstr.str();
          call_back.print_xmgrace( str );
        }
      }
  };

  // Gets an average of individuals accumulated over this generation
  // in kspace
  template< class t_Call_Back, class t_Object = typename t_Call_Back :: t_Object >
  class kPopAverage : public eoStatBase<t_Object>
  {
    protected:
      typedef typename std::complex< typename t_Object :: t_Type > t_complex;
      t_Call_Back &call_back;
      std::vector< t_complex > average;

    public:
      kPopAverage   ( t_Call_Back &_call_back )
                 : call_back( _call_back ) {}
      virtual ~kPopAverage() {}
      virtual std::string className(void) const { return "LaDa::kPopAverage"; }
      virtual void operator()( const eoPop<t_Object> &_pop )
      {
        t_unsigned pSize = _pop.size();
        typename eoPop<t_Object> :: const_iterator i_indiv = _pop.begin();
        typename eoPop<t_Object> :: const_iterator i_end = _pop.end();
        std::vector< t_complex > fourrier;
        
        // first iteration creates average
        call_back.fourrier_transform( *i_indiv, average );

        // then loops
        for( ++i_indiv; i_indiv != i_end; ++i_indiv )
        {
          call_back.fourrier_transform( *i_indiv, fourrier );
          std::transform( average.begin(), average.end(), fourrier.begin(),
                          average.begin(), std::plus<t_complex>() );
        }

        // prints stuff out
        {
          std::ostringstream sstr; 
          sstr << "PopAverage kspace: " << std::setw(5) << std::setprecision(2);
          typename std::vector<t_complex> :: iterator i_var = average.begin();
          typename std::vector<t_complex> :: iterator i_var_end = average.end();
          for( ; i_var != i_var_end; ++i_var )
            sstr << (*i_var / (t_real) pSize ) << " ";
          std::string str = sstr.str();
          call_back.print_xmgrace( str );
        }
      }
  };
  
  
  // Diversity measurment
  // \sum_{i in pop}\sum_{j > i} (\sigma_i - \sigma_j) / Npop / Npop-1 * 2 / Nspins
  template< class t_Call_Back, class t_Object = typename t_Call_Back :: t_Object >
  class Diversity : public eoStatBase<t_Object>
  {
    protected:
      t_Call_Back &call_back;
      bool keep_out_twins;

    public:
      Diversity   ( t_Call_Back &_call_back, bool _keepout = true ) 
                : call_back( _call_back ), keep_out_twins( _keepout ) {}
      virtual ~Diversity() {}
      virtual std::string className(void) const { return "LaDa::Diversity"; }
      virtual void operator()( const eoPop<t_Object> &_pop )
      {
        typename eoPop<t_Object> :: const_iterator i_indiv1 = _pop.begin();
        typename eoPop<t_Object> :: const_iterator i_end = _pop.end();
        typename eoPop<t_Object> :: const_iterator i_indiv2;
        typename t_Object :: const_iterator i_var1_begin, i_var1_end, i_var1, i_var2;
        t_real result = 0;
        t_unsigned N = 0;
        
        for( ; i_indiv1 != i_end; ++i_indiv1 )
        {
          i_var1_begin = i_indiv1->begin();
          i_var1_end = i_indiv1->end();
          for ( i_indiv2 = i_indiv1 + 1; i_indiv2 != i_end; ++i_indiv2 )
            if ( not ( *i_indiv1 == *i_indiv2 ) or not keep_out_twins )
            {
              ++N; 
              i_var1 = i_var1_begin; i_var2 = i_indiv2->begin();
              t_real d = 0;
              for (; i_var1 != i_var1_end; ++i_var1, ++i_var2 )
                d += ( *i_var1 - *i_var2 ) *  ( *i_var1 - *i_var2 );
              result += std::sqrt(d);
            }
        }

        N *= _pop.begin()->size();
        result /= static_cast<t_real>(N);

        // prints stuff out
        std::ostringstream sstr; 
        sstr << "Diversity: " << std::setw(10) << std::setprecision(5) << result;
        std::string str = sstr.str();
        call_back.print_xmgrace( str );
      }
  };

  // Census, discarding identical individuals
  template< class t_Call_Back, class t_Object = typename t_Call_Back :: t_Object >
  class TrueCensus : public eoStatBase<t_Object>
  {
    protected:
      t_Call_Back &call_back;

    public:
      TrueCensus   ( t_Call_Back &_call_back, bool _keepout = true ) 
                  : call_back( _call_back ) {}
      virtual ~TrueCensus() {}
      virtual std::string className(void) const { return "LaDa::TrueCensus"; }
      virtual void operator()( const eoPop<t_Object> &_pop )
      {
        typename eoPop<t_Object> :: const_iterator i_begin = _pop.begin();
        typename eoPop<t_Object> :: const_iterator i_indiv1 = i_begin;
        typename eoPop<t_Object> :: const_iterator i_end = _pop.end();
        typename eoPop<t_Object> :: const_iterator i_indiv2;
        t_unsigned N = 0;
        
        for( ; i_indiv1 != i_end; ++i_indiv1 )
        {
          for ( i_indiv2 = i_begin; i_indiv2 != i_end; ++i_indiv2 )
            if ( i_indiv1 != i_indiv2 and *i_indiv1 == *i_indiv2 )
              break;
          if ( i_indiv2 == i_end )
            ++N;
        }

        // prints stuff out
        std::ostringstream sstr; 
        sstr << "TruePopSize: " << std::setw(10) << std::setprecision(3)
             << N << " / " << _pop.size(); 
        std::string str = sstr.str();
        call_back.print_xmgrace( str );
      }
  };
} // namespace LaDa

#endif // _STATISTICS_H_
