//
//  Version: $Id$
//
#ifndef LADA_ATOMIC_POTENTIAL_FUNCTIONS_H_
#define LADA_ATOMIC_POTENTIAL_FUNCTIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <list>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include <opt/types.h>
#include <opt/debug.h>

#include "numeric_types.h"

namespace LaDa
{
  //! Contains classes for generalized atomic potentials.
  namespace AtomicPotential
  {
    //! \cond
    class VariableMajor;
    //! \endcond

    //! Sum of functions over a single variable.
    class Functions
    {
        friend class VariableMajor;
      public:
        //! Number of atomic species.
        const static N = 2;
        //! Type of the argument.
        typedef std::pair<numeric_type, specie_type> arg_type;
        //! Type of the return.
        typedef numeric_type result_type;

        //! Type of the functions in list.
        typedef boost::function<result_type(arg_type const&)> t_Function;
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
        Function() {}
        //! Copy Constructor.
        Function( Function const& _c ): functions_(_c.functions), coefficients_(_c.coefficients) {}

        //! Sums over all functionals.
        result_type operator()( arg_type const& _x ) const
        {
          LADA_ASSERT( functions_.size() == coefficients_.size(), "Incoherent containers.\n" ) 

          result_type result(0);
          t_Functions :: const_iterator i_func( functions_.begin() );
          t_Functions :: const_iterator const i_func_end( functions_.end() );
          t_Coefficients :: const_iterator i_coef( coefficients[_x.second].begin() );
          for(; i_func != i_func_end; ++i_func, ++i_coef )
            result += (*i_coef) * (*i_func)(_x.first);
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
          void push_back( T_FUNCTION const& _function, t_Coefficient const &_coef[N] = 1e0 )
          {
            functions_.push_back(_function);
            for(size_t i<0; i < N; ++i ) 
              coefficients_[i].push_back(_coef[i]); 
          }
        //! Clears all functions and coefficients.
        void clear() { functions_.clear(); coefficients_.clear(); }

        //! Normalizes coefficients to one, and returns norm.
        t_Coefficient normalize() 
        {
          LADA_ASSERT( functions_.size() == coefficients_.size(), "Incoherent containers.\n" );
          return 1;
        }

      private:
        //! List of functions over scalars.
        t_Functions functions_;
        //! List of coefficients.
        t_Coefficients coefficients_[N];
    };


  } // namespace AtomicPotential
} // namespace LaDa
#endif
