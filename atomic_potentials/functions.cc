//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef LADA_SEPN
# define LADA_SEPN 2
#endif

#include <boost/lambda/lambda.hpp>
#include "functions.h"

namespace LaDa
{
  namespace atomic_potential
  {

    Functions::result_type Functions::operator()( arg_type const& _x ) const
    {
#     ifdef LADA_DEBUG
        LADA_ASSERT( Functions::N * functions_.size() == coefficients_.size(),
                     "Incoherent containers.\n" ); 
#     endif

      result_type result(0);
      t_Functions :: const_iterator i_func( functions_.begin() );
      t_Functions :: const_iterator const i_func_end( functions_.end() );
      t_Coefficients :: const_iterator i_coef( coefficients_.begin() );
      for(; i_func != i_func_end; ++i_func, i_coef += Functions::N )
        result += (*(i_coef+_x.second)) * (*i_func)(_x.first);
      return result;
    }

    Functions::t_Coefficient Functions::normalize() 
    {
      namespace bl = boost::lambda; 
      LADA_ASSERT( Functions::N * functions_.size() == coefficients_.size(), 
                   "Incoherent containers.\n" );
      t_Coefficients::iterator i_first = coefficients_.begin();
      t_Coefficients::iterator const i_end = coefficients_.end();
      numeric_type const norm
      (
        std::sqrt
        ( 
          std::accumulate( i_first, i_end, numeric_type(0), bl::_1 + bl::_2 * bl::_2 )
        )
      );
      numeric_type const inv_norm( numeric_type(1) / norm );
      std::for_each(coefficients_.begin(), i_end, bl::_1 *= bl::constant(inv_norm));
      return norm;
    }

    std::ostream& operator<<( std::ostream &_stream,
                              Functions::const_iterator::reference const &_func )
    {
      _stream <<  "(" << _func[0];
      for(size_t i(1); i < Functions::N; ++i) _stream << ", " << _func[i];
      _stream << ")*" << _func.name() << " ";
      return _stream;
    }

    std::ostream& operator<<( std::ostream &_stream,
                              Functions::iterator::reference const &_func )
    {
      _stream <<  "(" << _func[0];
      for(size_t i(1); i < Functions::N; ++i) _stream << ", " << _func[i];
      _stream << ")*" << _func.name() << " ";
      return _stream;
    }

    std::ostream& operator<<( std::ostream &_stream, Functions const &_func )
    {
      Functions::const_iterator i_func( _func.begin() );
      Functions::const_iterator const i_func_end( _func.end() );
      _stream << *i_func; 
      for(++i_func; i_func != i_func_end; ++i_func) _stream << " + " << *i_func;

      return _stream;
    }

  } // namespace atomic_potential
} // namespace LaDa
