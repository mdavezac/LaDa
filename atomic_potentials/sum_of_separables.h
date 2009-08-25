//
//  Version: $Id$
//
#ifndef LADA_ATOMIC_POTENTIAL_SUM_OF_SEPARABLES_H_
#define LADA_ATOMIC_POTENTIAL_SUM_OF_SEPARABLES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/iterator/zip_iterator.hpp>

#include "separable.h"

namespace LaDa
{
  namespace atomic_potential
  {
    // Forward declaration.
    //! \cond
    namespace collapse
    {
      class CoordRange;
    }
    //! \endcond

    //! A separable function.
    class SumOfSeparables
    {
        friend class collapse::CoordRange;
      public:
        //! Argument type.
        typedef Separable::arg_type arg_type;
        //! Type of the return.
        typedef Separable::result_type result_type;

        //! Type of the functions in list.
        typedef Separable t_Function;
        //! Type of the function-list.
        typedef std::list<t_Function> t_Functions;
        //! Type of the list of coefficients.
        typedef result_type t_Coefficient;
        //! Type of the list of coefficients.
        typedef std::vector<t_Coefficient> t_Coefficients;

        //! Type of the iterator over functions and coefficients.
        typedef boost::zip_iterator
                <
                  boost::tuple<t_Functions::iterator, t_Coefficients::iterator>
                > iterator;
        //! Type of the iterator over functions and coefficients.
        typedef boost::zip_iterator
                <
                  boost::tuple<t_Functions::const_iterator, t_Coefficients::const_iterator>
                > const_iterator;
#       ifdef LADA_WITH_CONST
#         error LADA_WITH_CONST already defined
#       endif
        //! Iterates over coordinates first, then rank.
        class const_coord_range;
#       define LADA_WITH_CONST
#       include "sum_of_separables.coord_range.h"
        //! Iterates over coordinates first, then rank.
        class coord_range;
#       include "sum_of_separables.coord_range.h"

        
        //! Constructor.
        SumOfSeparables() {}
        //! Copy Constructor.
        SumOfSeparables   ( SumOfSeparables const& _c )
                        : functions_(_c.functions_), coefficients_(_c.coefficients_) {}

        //! Sums over all functionals.
        template<class T_CONTAINER> 
          result_type operator()( T_CONTAINER const& _x ) const
          {
            LADA_ASSERT( functions_.size() == coefficients_.size(), "Incoherent containers.\n" ) 
  
            result_type result(1);
            t_Functions :: const_iterator i_func( functions_.begin() );
            t_Functions :: const_iterator const i_func_end( functions_.end() );
            t_Coefficients :: const_iterator i_coef( coefficients_.begin() );
            for(; i_func != i_func_end; ++i_coef )
              result += (*i_coef) * (*i_func)(_x);
            return result;
          }


        //! Iterates over coordinates first, then rank.
        const_coord_range range() const { return const_coord_range(*this); }
        //! Iterates over coordinates first, then rank.
        coord_range range() { return coord_range(*this); }

        //! Returns iterator to functions and coefficients.
        iterator begin()
          { return iterator( boost::make_tuple(functions_.begin(), coefficients_.begin())); }
        //! Returns iterator to functions and coefficients.
        const_iterator begin() const
          { return const_iterator( boost::make_tuple(functions_.begin(), coefficients_.begin())); }
        //! Returns iterator to functions and coefficients.
        iterator end()
          { return iterator( boost::make_tuple(functions_.end(), coefficients_.end())); }
        //! Returns iterator to functions and coefficients.
        const_iterator end() const
          { return const_iterator( boost::make_tuple(functions_.end(), coefficients_.end())); }
        //! pushes a function and coefficient back.
        void push_back( Separable const& _function, t_Coefficient const &_coef = 1e0 )
          { functions_.push_back(_function); coefficients_.push_back(_coef); }
        //! Clears all functions and coefficients.
        void clear() { functions_.clear(); coefficients_.clear(); }
        //! Returns rank.
        size_t size() const { return coefficients_.size(); }
        //! Return number of coordinates.
        size_t nb_coordinates() const
        {
          size_t result(0);
          t_Functions::const_iterator i_sep( functions_.begin() );
          t_Functions::const_iterator const i_sep_end( functions_.end() );
          for(; i_sep != i_sep_end; ++i_sep)
            result = std::max( result, i_sep->size() );
        };

        //! Normalizes all separable functionals.
        void normalize()
        {
          LADA_ASSERT( functions_.size() == coefficients_.size(), "Incoherent containers.\n" ) 

          t_Functions :: iterator i_func( functions_.begin() );
          t_Functions :: iterator const i_func_end( functions_.end() );
          t_Coefficients :: iterator i_coef( coefficients_.begin() );
          for(; i_func != i_func_end; ++i_func, ++i_coef) 
            *i_coef *= i_func->normalize();
        }


      private:
        //! List of functions over scalars.
        t_Functions functions_;
        //! List of coefficients.
        t_Coefficients coefficients_;
    };


  } // namespace atomic_potential
} // namespace LaDa
#endif
