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
    return _stream << "    mse error: " << _b.get<0>() / _b.norm
                   << "    average error: " << _b.get<1>() / _b.norm
                   << "    maximum error: " << _b.get<2>();
  }
  void operator+=( ErrorTuple &_a, const ErrorTuple &_b )
  {
    _a.get<0>() += _b.get<0>();
    _a.get<1>() += _b.get<1>();
    _a.get<2>() = std::max( _a.get<2>(), _b.get<2>() );
    _a.norm += _b.norm;
  }
  std::ostream& operator<<( std::ostream &_stream, const NErrorTuple &_b )
  {
    return _stream << "    mse error: " << _b.get<0>() / _b.norm
                     << " ( " << 1e2 * _b.get<0>() /  _b.variance / _b.norm << "% )"
                   << "    average error: " << _b.get<1>() / _b.norm
                     << " ( " << 1e2 * _b.get<1>() / std::abs( _b.mean ) / _b.norm
                     << "% )"
                   << "    maximum error: " << _b.get<2>() 
                     << " ( " << 1e2 * _b.get<2>() / std::abs( _b.mean ) << "% )";
  }
}
