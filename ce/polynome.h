#ifndef _Polynome_H_
#define _Polynome_H_

#include "LaDaConfig.h"

#include <list>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>

#include <tinyxml/tinyxml.h>

#include <opt/function_base.h>
#include <opt/polynome.h>
#include <misc/types.h>

namespace LaDa
{
  namespace CE {

    //! \brief Adds a bit of XML and accounting to function::Polynome.
    //! \details The accounting comes from threre new static member variables,
    //!          Polynome::nb_eval, Polynome::nb_eval_grad, and Polynome::nb_eval_with_grad,
    //!          which are incremented by any instance of Polynome which makes a
    //!          call to Polynome::evaluate(), Polynome::evaluate_gradient(), and
    //!          Polynome::evaluate_with_gradient(), respectively.
    //!          Furthermore, this instance some XML saving/loading capability
    //!          for polynomials and monomials. 
    class Polynome : public function::Polynome<types::t_real>
    {
        //! Type of the base class
        typedef function::Polynome<types::t_real> t_Base;
        //! Type of the monomial components
        typedef function::Monome<types::t_real> t_Monome;
      public:
        //! Counts the calls to Polynome::evaluate().
        static types::t_unsigned nb_eval;
        //! Counts the calls to Polynome::evaluate_gradient().
        static types::t_unsigned nb_eval_grad;
        //! Counts the calls to Polynome::evaluate_with_gradient().
        static types::t_unsigned nb_eval_with_grad;

      public: 
        //! Constructor
        Polynome() : t_Base() {};
        //! Constructor and Initializer
        Polynome( types::t_int _nb ) : t_Base(_nb) {};
        //! Destructor
        virtual ~Polynome() {};
        

        //! Evaluates the polynome and increments Polynome::nb_eval.
        virtual types::t_real evaluate() 
          { ++nb_eval; return t_Base::evaluate();}
        //! Evaluates the gradient of the polynome and increments Polynome::nb_eval_grad.
        virtual void evaluate_gradient(types::t_real* const _grad) 
          { ++nb_eval_grad; return t_Base::evaluate_gradient( _grad ); }
        //! Evaluates the polynome and its gradient and increments Polynome::nb_eval_with_grad.
        virtual types::t_real evaluate_with_gradient(types::t_real* const _grad) 
          { ++nb_eval_with_grad; return t_Base::evaluate_with_gradient( _grad ); }

        //! \brief Resizes the number of variables according to the largest
        //!        monomial variable index.
        void resize_variables();
        //! Dumps a polynomial to an XML node.
        void print_xml( TiXmlElement* const node ) const;
        //! Dumps a monomial to an XML node.
        void print_xml( const t_Monome &_m, TiXmlElement* const node ) const;
        //! Load a monomials from an XML node.
        bool Load( t_Monome &_m, TiXmlElement* const node );
        //! Load a polynomial from an XML node.
        bool Load( TiXmlElement* const node );

        #ifdef LADA_DEBUG
          //! Doesn't do anything...
          bool is_consistent() const;
        #endif 
    };

    
  } // namespace CE
} // namespace LaDa
#endif
