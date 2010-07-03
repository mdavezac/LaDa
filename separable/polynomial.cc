#include "LaDaConfig.h"

#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include "polynomial.h"

namespace LaDa
{
  namespace Separable
  {

    Polynomial :: set_range( types::t_int _min, types::t_int _max )
    {
      if( _max < _min ) std::swap( _min, _max );
      if( _max == _min ) ++_max;
      using namespace boost::lambda;
      min = _min; 
      max = _max;
      container.resize( _max - _min );
      types::t_int i = _min-1;
      std::generate( container.begin(), container.end(),
                     bind( &Polynomial::degree, _1 ) = ( ++var(i) ) );
    }

  } // end of Separable namespace
} // namespace LaDa
