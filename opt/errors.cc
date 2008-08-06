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
    return _stream << "    mse error: " << _b.variance()
                   << "    average error: " << _b.mean()
                   << "    maximum error: " << _b.max();
  }
  ErrorTuple& ErrorTuple::operator+=( const ErrorTuple &_b )
  {
    get<0>() += _b.get<0>();
    get<1>() += _b.get<1>();
    get<2>() = std::max( get<2>(), _b.get<2>() );
    norm += _b.norm;
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
}
