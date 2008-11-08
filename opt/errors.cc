//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/lambda/bind.hpp>
#include "errors.h"

namespace LaDa
{
  namespace opt
  {
    std::ostream& operator<<( std::ostream &_stream, const ErrorTuple &_b )
    {
      return _stream << "    mse error: " << _b.variance()
                     << "    average error: " << _b.mean()
                     << "    maximum error: " << _b.max();
    }
    ErrorTuple& ErrorTuple::operator+=( const ErrorTuple &_b )
    {
      get<0>() += _b.get<0>();
      get<1>() += _b.get<1>();
      get<2>() = std::max( get<2>(), _b.get<2>() );
      norm_ += _b.norm_;
      return *this;
    }
    std::ostream& operator<<( std::ostream &_stream, const NErrorTuple &_b )
    {
      return _stream << "    mse error: " << _b.ErrorTuple::variance()
                       << " ( " << 1e2 * _b.variance() << "% )"
                     << "    average error: " << _b.ErrorTuple::mean()
                       << " ( " << 1e2 * _b.mean() << "% )"
                     << "    maximum error: " << _b.ErrorTuple::max()
                       << " ( " << 1e2 * _b.max() << "% )";
    }
    NErrorTuple mean_n_var( const std::vector<Crystal::Structure> &_strs )
    {
      namespace bl = boost::lambda;
      std::vector<types::t_real> t, w;
      std::transform
      (
        _strs.begin(), _strs.end(), std::back_inserter( t ),
        bl::bind( &Crystal::Structure::energy, bl::_1 )
      );
      std::transform
      (
        _strs.begin(), _strs.end(), std::back_inserter( w ),
        bl::bind( &Crystal::Structure::weight, bl::_1 )
      );
      return mean_n_var( t, w );
    }
  }
} // namespace LaDa
