//
// Monome class: a single term made up of a coefficient and an array
// of variables. The latters are integers which are meant to refer to
// the variable array of a Polynome object
//
// Polynome class: contains an array of Monome objects, as well as an
// array of variables which make up the monomes 
//
// Constrained minimization can be achieved via the OPT++ package
// and the opt_minimize.h template class
// see lamarck.cc and functional_builder.cc for an example
//                 
#ifndef _Polynome_H_
#define _Polynome_H_

#include <list>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <tinyxml/tinyxml.h>
#include <opt/opt_function_base.h>
#include <opt/opt_polynome.h>
#include <opt/types.h>

namespace VA_CE {

  class Polynome : public function::Polynome<types::t_real>
  {
    public:
      static types::t_unsigned nb_eval, nb_eval_grad, nb_eval_with_grad;

    public: 
      Polynome() : function::Polynome<types::t_real>() {};
      Polynome( types::t_int _nb ) : function::Polynome<types::t_real>(_nb) {};
      virtual ~Polynome() {};
      

      // evaluations
      virtual types::t_real evaluate() 
        { ++nb_eval; return function::Polynome<types::t_real>::evaluate();}
      virtual void evaluate_gradient(types::t_real* const _grad) 
        { ++nb_eval_grad; return function::Polynome<types::t_real>::evaluate_gradient( _grad ); }
      virtual types::t_real evaluate_with_gradient(types::t_real* const _grad) 
        { ++nb_eval_with_grad; return function::Polynome<types::t_real>::evaluate_with_gradient( _grad ); }

      void resize_variables();
      void print_xml( TiXmlElement* const node ) const;
      void print_xml( const function::Monome<> &_m, TiXmlElement* const node ) const;
      bool Load( function::Monome<> &_m, TiXmlElement* const node );
      bool Load( TiXmlElement* const node );

      #ifdef _DEBUG_LADA_ 
        bool is_consistent() const;
      #endif // _DEBUG_LADA_ 
  };

  
} // namespace VA_CE
#endif
