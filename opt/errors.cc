//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "errors.h"

namespace opt
{
  std::ostream& operator<<( std::ostream &_stream, const ErrorTuple &_b )
  {
    return _stream << "    mse error: " << _b.get<0>()
                   << "    average error: " << _b.get<1>()
                   << "    maximum error: " << _b.get<2>();
  }
  void operator+=( ErrorTuple &_a, const ErrorTuple &_b )
  {
    _a.get<0>() += _b.get<0>();
    _a.get<1>() += _b.get<1>();
    _a.get<2>() = std::max( _a.get<2>(), _b.get<2>() );
  }
  void operator/=( ErrorTuple &_a, const types::t_unsigned _N )
  {
    _a.get<0>() /= types::t_real( N );
    _a.get<1>() /= types::t_real( N );
  }
  std::ostream& operator<<( std::ostream &_stream, const NErrorTuple &_b )
  {
    return _stream << "    mse error: " << _b.get<0>()
                     << " ( " << 1e2 * _b.get<0>() /  _b.variance << "% )"
                   << "    average error: " << _b.get<1>()
                     << " ( " << 1e2 * _b.get<1>() / std::abs( _b.mean ) << "% )"
                   << "    maximum error: " << _b.get<2>()
                     << " ( " << 1e2 * _b.get<2>() / std::abs( _b.mean ) << "% )";
  }
}
