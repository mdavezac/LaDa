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

#include <opt/types.h>
#include <opt/debug.h>

namespace LaDa
{
  namespace Gsl
  {
    //! \brief Wraps gsl_vector around a container.
    //! \warning Expects contiguous memory. 
    class Vector
    {
      protected:
        //! Defines the memory block.
        gsl_block bloc;
        //! Defines the gsl vector.
        gsl_vector vector;

      public:
        //! Constructor and Initializer.
        template< class T_CONTAINER > Vector( T_CONTAINER& _vec );
        //! Destructor.
        virtual ~Vector() {};

        //! Returns a pointer to the gsl_vector.
        operator gsl_vector*() { return &vector; }
        //! Returns a const pointer to the gsl_vector.
        operator const gsl_vector* () const { return &vector; }
    };
    //! \brief Wraps gsl_maxtrix around a container.
    //! \warning Expects contiguous memory. See gsl_matrix documentation for
    //!          layout details.
    class Matrix
    {
      protected:
        //! Defines the memory block.
        gsl_block bloc;
        //! Defines the gsl maxtrix.
        gsl_matrix matrix;

      public:
        //! Constructor and Initializer.
        Matrix( types::t_int _nrow, std::vector<types::t_real>& _vec );
        //! Destructor.
        ~Matrix() {};

        //! Returns a pointer to the gsl_matrix.
        operator gsl_matrix*() { return &matrix; }
        //! Returns a const pointer to the gsl_matrix.
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
