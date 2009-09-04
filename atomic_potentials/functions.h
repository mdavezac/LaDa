//
//  Version: $Id$
//
#ifndef LADA_ATOMIC_POTENTIAL_FUNCTIONS_H_
#define LADA_ATOMIC_POTENTIAL_FUNCTIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <list>
#include <iostream>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <opt/types.h>
#include <opt/debug.h>

#include "numeric_types.h"

namespace LaDa
{
  //! Contains classes for generalized atomic potentials.
  namespace atomic_potential
  {
    class Functions
    {
      public:
        //! Number of atomic species.
        const static size_t N = 2;
        //! Type of the argument.
        typedef std::pair<numeric_type, specie_type> arg_type;
        //! Type of the return.
        typedef numeric_type result_type;

        //! Type of the functions in list.
        typedef boost::function<result_type(arg_type::first_type const&)> t_Function;
        //! Type of the function-list.
        typedef std::list<t_Function> t_Functions;
        //! Type of the list of coefficients.
        typedef result_type t_Coefficient;
        //! Type of the list of coefficients.
        typedef std::vector<t_Coefficient> t_Coefficients;

        //! Type of the iterator over functions and coefficients.
        class iterator;
        //! Type of the iterator over functions and coefficients.
        class const_iterator;

        // includes iterator definitions here.
#       if defined(LADA_WITH_CONST) or defined(LADA_ITERATOR_NAME)
#         error LADA_WITH_CONST and LADA_ITERATOR_NAME already defined.
#       endif
#       define LADA_ITERATOR_NAME iterator
#       define LADA_WITH_CONST 
#       include "functions.iterator.h"
#       define LADA_ITERATOR_NAME const_iterator
#       define LADA_WITH_CONST const
#       include "functions.iterator.h"

        //! Constructor.
        Functions() {}
        //! Copy Constructor.
        Functions   ( Functions const& _c )
                  : functions_(_c.functions_), coefficients_(_c.coefficients_) {}

        //! Sums over all functionals.
        result_type operator()( arg_type const& _x ) const
        {
#         ifdef LADA_DEBUG
            LADA_ASSERT( Functions::N * functions_.size() == coefficients_.size(),
                         "Incoherent containers.\n" ); 
#         endif

          result_type result(0);
          t_Functions :: const_iterator i_func( functions_.begin() );
          t_Functions :: const_iterator const i_func_end( functions_.end() );
          t_Coefficients :: const_iterator i_coef( coefficients_.begin() );
          for(; i_func != i_func_end; ++i_func, i_coef += Functions::N )
            result += (*(i_coef+_x.second)) * (*i_func)(_x.first);
          return result;
        }

        //! Returns iterator to functions and coefficients.
        iterator begin()
          { return iterator( functions_.begin(), coefficients_.begin()); }
        //! Returns iterator to functions and coefficients.
        const_iterator begin() const
          { return const_iterator( functions_.begin(), coefficients_.begin()); }
        //! Returns iterator to functions and coefficients.
        iterator end()
          { return iterator( functions_.end(), coefficients_.end()); }
        //! Returns iterator to functions and coefficients.
        const_iterator end() const
          { return const_iterator( functions_.end(), coefficients_.end()); }
        //! pushes a function and coefficient back.
        template< class T_FUNCTION >
          void push_back( T_FUNCTION const& _function, numeric_type const _coef[N] )
          {
            functions_.push_back(_function);
            for(size_t i(0); i < N; ++i ) 
              coefficients_.push_back(_coef[i]); 
          }
        //! Clears all functions and coefficients.
        void clear() { functions_.clear(); coefficients_.clear(); }

        //! Normalizes coefficients to one, and returns norm.
        t_Coefficient normalize() 
        {
          LADA_ASSERT( functions_.size() == coefficients_.size(), "Incoherent containers.\n" );
          return 1;
        }

        //! Returns the number of functions.
        size_t size() const { return functions_.size(); }

      private:
        //! List of functions over scalars.
        t_Functions functions_;
        //! List of coefficients.
        t_Coefficients coefficients_;
    };

    inline std::ostream& operator<<( std::ostream &_stream, Functions const &_func )
    {
      Functions::const_iterator i_func( _func.begin() );
      Functions::const_iterator const i_func_end( _func.end() );
      for(; i_func != i_func_end; ++i_func)
      {
        _stream << "(" << (*i_func)[0];
        for(size_t i(1); i < _func.N; ++i)
          _stream << ", " << (*i_func)[i];
        _stream << ") ";
      }

      return _stream;
    }


  } // namespace atomic_potential
} // namespace LaDa
#endif
