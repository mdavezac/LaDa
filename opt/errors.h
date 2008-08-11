//
//  Version: $Id$
//

#ifndef _OPT_ERRORS_H_
#define _OPT_ERRORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/tuple/tuple.hpp>

#include <vector>

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/algorithms.h>
#include <crystal/structure.h>

namespace opt
{
  //! An error tuple.
  class ErrorTuple : public boost::tuple<types::t_real, types::t_real, types::t_real> 
  {  
    protected:
      //! Base class.
      typedef boost::tuple<types::t_real, types::t_real, types::t_real> t_Base;
      //! Number of quantities over which to average, or sum of their weights in average.
      types::t_real norm;
    public:
      //! Default constructor.
      ErrorTuple() : t_Base(0,0,0), norm(0) {}
      //! \brief Constructor.
      //! \details ErrorTuple::norm is set to one.
      //! \param[in] _a is the variance.
      //! \param[in] _b is the mean.
      //! \param[in] |_c| is the max.
      ErrorTuple  ( types::t_real _a, types::t_real _b, types::t_real _c )
                : t_Base( _a, _b, _c ), norm(1)  {} 
      //! \brief Constructor. 
      //! \details The variance is set to \a _a * \a _a * \a _b, 
      //!          the mean to \a _a * \a _b, and the max to |\a _a|.
      //! \param[in] _a is the error.
      //! \param[in] _b is the weight.
      ErrorTuple  ( types::t_real _a, types::t_real _b = 1e0 )
                : t_Base( _a * _a * _b, std::abs(_a) * _b, std::abs(_a) ),
                  norm(_b){} 
      //! Copy constructor, 
      ErrorTuple  ( const ErrorTuple &_e ) : t_Base( _e ), norm(_e.norm) {} 
      //! returns variance.
      types::t_real variance() const { return get<0>() / norm; }
      //! returns mean.
      types::t_real mean() const { return get<1>() / norm; }
      //! returns max.
      types::t_real max() const { return get<2>(); }
      //! Sums errors.
      ErrorTuple& operator+=( const ErrorTuple &_b );
  };

  //! A pair of error tuple for prediction/training.
  typedef std::pair<ErrorTuple, ErrorTuple> t_ErrorPair;

  //! A normalized error tuple.
  struct NErrorTuple : public ErrorTuple
  {
    //! Base class.
    typedef ErrorTuple t_Base;
    //! Variance.
    types::t_real variance_;
    //! Mean.
    types::t_real mean_;
 
    public:
      //! Default constructor.
      NErrorTuple() : t_Base(), variance_(0), mean_(0) {}
      //! \brief Constructor. Normalization is set to zero.
      //! \details ErrorTuple::norm is set to one.
      //! \param[in] _a is the variance.
      //! \param[in] _b is the mean.
      //! \param[in] _c is the max.
      NErrorTuple  ( types::t_real _a, types::t_real _b, types::t_real _c )
                 : t_Base( _a, _b, _c ), variance_(0), mean_(0) {} 
      //! \brief Constructor. Normalization is set to zero.
      //! \details The variance is set to \a _a * \a _a * \a _b, 
      //!          the mean to \a _a * \a _b, and the max to \a _a.
      //! \param[in] _a is the error.
      //! \param[in] _b is the weight.
      NErrorTuple  ( types::t_real _a, types::t_real _b = 1e0 )
                 : t_Base( _a, _b ), variance_(0), mean_(0) {}
      //! Copy constructor.
      NErrorTuple  ( const NErrorTuple &_e )
                 : t_Base( _e ), variance_( _e.variance_ ), mean_( _e.mean_ ) {}
      //! Copies error tuple and sum of weights into normalized error tuple.
      const NErrorTuple& operator=( const ErrorTuple &_e )
       { *( (ErrorTuple*)this ) = _e; return *this; }
      //! returns variance.
      types::t_real variance() const { return ErrorTuple::variance() / variance_; }
      //! returns mean.
      types::t_real mean() const { return ErrorTuple::mean() / std::abs(mean_); }
      //! returns max.
      types::t_real max() const { return ErrorTuple::max() / std::abs(mean_); }
      //! return normalization coef of mean.
      types::t_real& nmean() { return mean_; }
      //! return normalization coef of mean.
      const types::t_real& nmean() const { return mean_; }
      //! return normalization coef of variance.
      types::t_real& nvariance() { return variance_; }
      //! return normalization coef of variance.
      const types::t_real& nvariance() const { return variance_; }
  };
  //! computes mean and variance of data
  template< class T_VECTOR >
  NErrorTuple mean_n_var( const T_VECTOR &_targets, const T_VECTOR &_weights );
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
    typename T_VECTOR::const_iterator i_weight = _weights.begin();
    __ASSERT( _targets.size() != _weights.size(),
              "Inconsistent vector sizes.\n" );
    for(; i_target != i_target_end; ++i_target, ++i_weight )
    {
      mean += (*i_target) * (*i_weight);
      norm += (*i_weight);
    }
    nerror.mean_ = types::t_real(mean) / types::t_real(norm);
    i_target = _targets.begin();
    i_weight = _weights.begin();
    for(; i_target != i_target_end; ++i_target, ++i_weight )
      square +=   ( (*i_target) - mean / norm )
                * ( (*i_target) - mean / norm )  
                * (*i_weight);
    nerror.variance_ = types::t_real(square) / types::t_real(norm);
    return nerror;
    __DEBUGTRYEND(, "Error in opt::mean_n_var().\n" )
  }

  // Computes mean and variance of energies in structure set.
  NErrorTuple mean_n_var( const std::vector<Crystal::Structure> &_strs );

}

#endif
