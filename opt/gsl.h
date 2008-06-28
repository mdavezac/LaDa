//
//  Version: $Id$
//
#ifndef _OPT_GSL_H_
#define _OPT_GSL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "types.h"
#include "debug.h"

namespace Gsl
{
  class Vector
  {
    protected:
      gsl_block bloc;
      gsl_vector vector;

    public:
      Vector( std::vector<types::t_real>& _vec );
      ~Vector() {};

      operator gsl_vector*() { return &vector; }
      operator const gsl_vector* const () const { return &vector; }
  };
  class Matrix
  {
    protected:
      gsl_block bloc;
      gsl_matrix matrix;

    public:
      Matrix( types::t_int _nrow, std::vector<types::t_real>& _vec );
      ~Matrix() {};

      operator gsl_matrix*() { return &matrix; }
      operator const gsl_matrix* const () const { return &matrix; }
  };
}
#endif
