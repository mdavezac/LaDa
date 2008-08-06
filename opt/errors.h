//
//  Version: $Id$
//

#ifndef _OPT_ERRORS_H_
#define _OPT_ERRORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/tuple/tuple.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/algorithms.h>

namespace opt
{
  //! An error tuple.
  class ErrorTuple : public boost::tuple<types::t_real, types::t_real, types::t_real> 
  {  
    //! Base class.
    typedef boost::tuple<types::t_real, types::t_real, types::t_real> t_Base;
    public:
      types::t_real norm;
      ErrorTuple() : t_Base(0,0,0), norm(0) {}
      ErrorTuple  ( types::t_real _a, types::t_real _b, types::t_real _c )
                : t_Base( _a, _b, _c ), norm(1)  {} 
      ErrorTuple  ( types::t_real _a, types::t_real _b )
                : t_Base( _a * _a * _b, _a * _b, std::max( _a, get<2>() ) ),
                  norm(_b){} 
      ErrorTuple  ( types::t_real _a )
                : t_Base( _a * _a, _a, std::max( _a, get<2>() ) ), norm(1) {} 
      ErrorTuple  ( const ErrorTuple &_e ) : t_Base( _e ), norm(_e.norm) {} 
  };
  //! A normalized error tuple.
  struct NErrorTuple : public ErrorTuple
  {
    //! Base class.
    typedef ErrorTuple t_Base;
    public:
      //! Variance.
      types::t_real variance;
      //! Mean.
      types::t_real mean;
 
      NErrorTuple() : t_Base(), variance(0), mean(0) {}
      NErrorTuple  ( types::t_real _a, types::t_real _b, types::t_real _c )
                 : t_Base( _a, _b, _c ), variance(0), mean(0) {} 
      NErrorTuple  ( types::t_real _a, types::t_real _b )
                 : t_Base( _a * _a * _b, _a * _b, std::max( _a, get<2>() ) ),
                   variance(0), mean(0) {}
      NErrorTuple  ( types::t_real _a )
                 : t_Base( _a * _a, _a, std::max( _a, get<2>() ) ),
                   variance(0), mean(0) {}
      NErrorTuple  ( const NErrorTuple &_e )
                 : t_Base( _e ), variance( _e.variance ), mean( _e.mean ) {}
      const NErrorTuple& operator=( const ErrorTuple &_e )
      { 
        get<0>() = _e.get<0>(); get<1>() = _e.get<1>(); get<2>() = _e.get<2>(); 
        return *this;
      }
  };
  //! computes mean and variance of data
  template< class T_VECTOR >
  NErrorTuple mean_n_var( const T_VECTOR &_targets, const T_VECTOR &_weights );
  //! Sums errors.
  void operator+=( ErrorTuple &_a, const ErrorTuple &_b );
  //! Outputs an error tuple.
  std::ostream& operator<<( std::ostream &_stream, const ErrorTuple &_b );
  //! Outputs a normalized error tuple.
  std::ostream& operator<<( std::ostream &_stream, const NErrorTuple &_b );

  template< class T_VECTOR >
  NErrorTuple mean_n_var( const T_VECTOR &_targets, const T_VECTOR &_weights )
  {
    __DEBUGTRYBEGIN
    typename T_VECTOR::value_type norm( 0 ), square(0), mean(0);
    NErrorTuple nerror;
    typename T_VECTOR::const_iterator i_target = _targets.begin();
    typename T_VECTOR::const_iterator i_target_end = _targets.end();
    typename T_VECTOR::const_iterator i_weight = _weights.end();
    __ASSERT( _targets.size() != _weights.size(),
              "Inconsistent vector sizes.\n" );
    for(; i_target != i_target_end; ++i_target )
    {
      mean += (*i_target) * (*i_weight);
      square += (*i_target) * (*i_target) * (*i_weight);
      norm += (*i_weight);
    }
    nerror.mean = types::t_real(mean) / types::t_real(norm);
    nerror.variance = types::t_real(square - mean*mean) / types::t_real(norm);
    return nerror;
    __DEBUGTRYEND(, "Error in opt::mean_n_var().\n" )
  }

}

#endif
