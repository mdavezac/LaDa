//
//  Version: $Id$
//
#ifndef _OPT_GSL_H_
#define _OPT_GSL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include<boost/type_traits/is_same.hpp>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "types.h"
#include "debug.h"

namespace LaDa
{
  namespace Gsl
  {
    class Vector
    {
      protected:
        gsl_block bloc;
        gsl_vector vector;

      public:
        template< class T_CONTAINER > Vector( T_CONTAINER& _vec );
        ~Vector() {};

        operator gsl_vector*() { return &vector; }
        operator const gsl_vector* () const { return &vector; }
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
        operator const gsl_matrix* () const { return &matrix; }
    };
    template< class T_CONTAINER >
      Vector::Vector( T_CONTAINER &_vec )
      {
        if( not boost::is_same< typename T_CONTAINER::value_type, double>::value )
         __THROW_ERROR("Types are not equivalent.\n" )
        vector.size = _vec.size();
        vector.stride = 1;
        vector.data = &_vec[0];
        vector.block = &bloc;
        vector.owner = 0;
        bloc.size = vector.size;
        bloc.data = vector.data;
      }
  }
} // namespace LaDa
#endif
