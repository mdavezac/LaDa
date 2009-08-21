//
//  Version: $Id$
//
#ifndef LADA_ATOMIC_POTENTIAL_SUM_OF_SEPARABLES_H_
#define LADA_ATOMIC_POTENTIAL_SUM_OF_SEPARABLES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "separable.h"

namespace LaDa
{
  namespace AtomicPotential
  {
    // Forward declaration.
    //! \cond
    namespace collapse
    {
      class VariableMajorRange;
    }
    //! \endcond

    //! A separable function.
    class SumOfSeparables
    {
      public:
        //! Argument type.
        typedef Separable::arg_type arg_type;
        //! Type of the return.
        typedef Separable::result_type result_type;

        //! Type of the functions in list.
        typedef Functions t_Function;
        //! Type of the function-list.
        typedef std::list<t_Function> t_Functions;
        //! Type of the list of coefficients.
        typedef result_type t_Coefficient;
        //! Type of the list of coefficients.
        typedef boost::numeric::ublas::vector<t_Coefficient> t_Coefficients;

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
        
        //! Constructor.
        SumOfSeparables() {}
        //! Copy Constructor.
        SumOfSeparables   ( SumOfSeparables const& _c )
                        : functions_(_c.functions), coefficients_(_c.coefficients) {}

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

        //! Returns iterator to functions and coefficients.
        iterator begin()
          { return iterator( boost::make_tuple(functions_.begin(), coefficients.begin())); }
        //! Returns iterator to functions and coefficients.
        const_iterator begin() const
          { return iterator( boost::make_tuple(functions_.begin(), coefficients.begin())); }
        //! Returns iterator to functions and coefficients.
        iterator end()
          { return iterator( boost::make_tuple(functions_.end(), coefficients.end())); }
        //! Returns iterator to functions and coefficients.
        const_iterator end() const
          { return iterator( boost::make_tuple(functions_.end(), coefficients.end())); }
        //! pushes a function and coefficient back.
        template< class T_FUNCTION >
          void push_back( T_FUNCTION const& _function, t_Coefficient const &_coef = 1e0 )
            { functions_.push_back(_function); coefficients_.push_back(_coef); }
        //! Clears all functions and coefficients.
        void clear() { functions_.clear(); coefficients_.clear(); }
        //! Returns rank.
        size_t size() const { return coefficients_.size(); }

        //! Normalizes all separable functionals.
        void normalize()
        {
          LADA_ASSERT( functions_.size() == coefficients_.size(), "Incoherent containers.\n" ) 

          t_Functions :: iterator i_func( functions_.begin() );
          t_Functions :: iterator const i_func_end( functions_.end() );
          t_Coefficients :: const_iterator i_coef( coefficients_.begin() );
          for(; i_func != i_func_end; ++i_func, ++i_coef) 
            *i_coef *= i_func->normalize();
        }


      private:
        //! List of functions over scalars.
        t_Functions functions_;
        //! List of coefficients.
        t_Coefficients coefficients_;
    };


  } // namespace AtomicPotential
} // namespace LaDa
#endif
