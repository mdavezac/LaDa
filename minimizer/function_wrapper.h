#ifndef _LADA_MINIMIZER_FUNCTION_WRAPPER_H_
#define _LADA_MINIMIZER_FUNCTION_WRAPPER_H_

#include "LaDaConfig.h"

#include<algorithm>

namespace LaDa
{
  namespace function
  {

    //! Wraps around a function::Base for use with minimizers.
    template<class T_FUNCTION >
    class MinimizerWrapper
    {
      public:
        //! Type of the functional.
        typedef T_FUNCTION t_Function;
        //! Return value.
        typedef typename T_FUNCTION :: t_Type t_Return; 
        //! Argument
        typedef typename T_FUNCTION :: t_Container t_Arg;

        //! \brief Constructor and Initializer
        //! \details can assign an external container as Base::variables
        MinimizerWrapper( t_Function& _function ) : function_( &_function ) {}
        //! Copy Constructor
        MinimizerWrapper( const MinimizerWrapper &_c ) : function_( _c.function_ ) {}
        
        //! Returns arguments
        t_Arg& argument() { return *(function_->variables); }
        //! Returns arguments
        const t_Arg& argument() const { return *function_->variables; }

        //! Evaluates.
        t_Return operator()( t_Return* const _arg ) const
        { 
          std::copy( _arg, _arg + function_->size(), function_->begin() );
          function_->evaluate();
        }
        //! Gradient value;
        void gradient( t_Return* const _arg, t_Return *_grad ) const
        {
          std::copy( _arg, _arg + function_->size(), function_->begin() );
          function_->evaluate_gradient( _grad ); 
        }

      protected:
        //! The function to wrap.
        t_Function* function_;
    };

    template< class T_FUNCTION >
      MinimizerWrapper<T_FUNCTION> minimizer_wrapper( T_FUNCTION& _function )
        { return MinimizerWrapper<T_FUNCTION>( &_function ); }

  }
} // namespace LaDa

#endif // _OPT_FUNCTION_BASE_H_
