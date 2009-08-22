//
//  Version: $Id: collapse.h 1270 2009-08-17 03:11:58Z davezac $
//
#ifndef LADA_ATOMIC_POTENTIAL_NUMERIC_TYPES_H_
#define LADA_ATOMIC_POTENTIAL_NUMERIC_TYPES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <opt/types.h>

namespace LaDa
{
  namespace atomic_potential
  {
    //! Basic numeric type.
    typedef types::t_real numeric_type;
    //! Basic specie type.
    typedef size_t specie_type;
    //! Numeric vector type.
    typedef boost::numeric::ublas::vector<numeric_type> vector_type;
    //! Numeric matrix type.
    typedef boost::numeric::ublas::matrix<numeric_type> matrix_type;
  } // namespace atomic_potential
} // namespace LaDa
#endif
