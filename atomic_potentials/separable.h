//
//  Version: $Id$
//
#ifndef LADA_ATOMIC_POTENTIAL_SEPARABLE_H_
#define LADA_ATOMIC_POTENTIAL_SEPARABLE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "functions.h"

namespace LaDa
{
  namespace atomic_potential
  {
    //! A separable function.
    class Separable
    {
        friend class VariableMajor;
      public:
        //! Argument type.
        typedef std::vector<Functions::arg_type> arg_type;
        //! Type of the return.
        typedef Functions::result_type result_type;

        //! Type of the functions in list.
        typedef Functions t_Function;
        //! Type of the function-list.
        typedef std::list<t_Function> t_Functions;

        //! Type of the iterator over functions and coefficients.
        typedef t_Functions::iterator iterator;
        //! Type of the iterator over functions and coefficients.
        typedef t_Functions::const_iterator const_iterator;

        //! Constructor.
        Separable() {}
        //! Copy Constructor.
        Separable   ( Separable const& _c ) : functions_(_c.functions_) {}

        //! Sums over all functionals.
        template<class T_CONTAINER> 
          result_type operator()( T_CONTAINER const& _x ) const
          {
            LADA_ASSERT( functions_.size() <= _x.size(), "Incoherent containers.\n" ) 
  
            result_type result(1);
            t_Functions :: const_iterator i_func( functions_.begin() );
            t_Functions :: const_iterator const i_func_end( functions_.end() );
            typename T_CONTAINER :: const_iterator i_x( _x.begin() );
            for(; i_func != i_func_end; ++i_func, ++i_x)
              result *= (*i_func)(*i_x);
            return result;
          }

        //! Returns iterator to functions and coefficients.
        iterator begin() { return functions_.begin(); }
        //! Returns iterator to functions and coefficients.
        const_iterator begin() const { return functions_.begin(); }
        //! Returns iterator to functions and coefficients.
        iterator end() { return functions_.end(); }
        //! Returns iterator to functions and coefficients.
        const_iterator end() const { return functions_.end(); }
        //! pushes a function and coefficient back.
        void push_back( t_Function const& _function )
          { functions_.push_back(_function); }
        //! Clears all functions and coefficients.
        void clear() { functions_.clear(); }
        //! Returns number of variables.
        size_t size() const { return functions_.size(); }
        //! Returns number of variables.
        size_t nb_coordinates() const { return size(); }


        result_type normalize()
        {
          result_type result(1);
          t_Functions :: iterator i_func( functions_.begin() );
          t_Functions :: iterator const i_func_end( functions_.end() );
          for(; i_func != i_func_end; ++i_func ) result *= i_func->normalize();
          return result;
        }

      private:
        //! List of functions over scalars.
        t_Functions functions_;
    };


  } // namespace atomic_potential
} // namespace LaDa
#endif
