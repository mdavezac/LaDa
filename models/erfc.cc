#include "LaDaConfig.h"

#include <boost/math/special_functions/erf.hpp>
extern "C" double FC_FUNC_(boost_erfc, BOOST_ERFC)( const double *const _in )
  { return boost::math::erfc( *_in ); }
