//
// Monome class: a single term made up of a coefficient and an array
// of variables. The latters are usually integers which are meant to refer to
// the variable array of a Polynome object
//
//                 
#ifndef _OPT_MONOME_H_
#define _OPT_MONOME_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <list>
#include <algorithm>
#include <numeric>
#include <iostream>

#include "opt/types.h"

namespace function {

  template<class _TYPE, class _TERM> class Polynome;
  template<class _COEF=types::t_real, class _TERM=types::t_int >
  class Monome
  {
    friend class Polynome<_COEF, _TERM>;
    public:
      typedef _COEF COEF_TYPE;
      typedef _TERM TERM_TYPE;
      typedef typename std::list<_TERM> TERM_CONTAINER_TYPE;
      typedef typename std::list<_TERM> :: iterator iterator;
      typedef typename std::list<_TERM> :: const_iterator const_iterator;

    protected:
      COEF_TYPE coefficient;
      std::list<TERM_TYPE> terms;

    public:
      Monome() : coefficient(0) {};
      Monome(const COEF_TYPE &_coef) : coefficient (_coef) {};
      Monome(const Monome<_COEF, _TERM> &_m) : coefficient (_m.coefficient), terms(_m.terms) {};
      ~Monome(){};

      // initialisation
      void set_coefficient( const COEF_TYPE &_c ) 
        { coefficient = _c; };
      COEF_TYPE get_coefficient() const
        { return coefficient; };
      size_t order() const
        { return terms.empty() ? 0 : terms.size(); }
      typename std::list<_TERM> :: const_iterator terms_begin() const
        { return terms.begin(); }
      typename std::list<_TERM> :: iterator terms_begin() 
        { return terms.begin(); }
      typename std::list<_TERM> :: const_iterator terms_end() const
        { return terms.end(); }
      typename std::list<_TERM> :: iterator terms_end()
        { return terms.end(); }

        // adds a term to to the monome. 
        // this is an index to the variables defined in a Polynome object
        // first add_term checks for memory availability
        // then it places the value of _term in a sorted list
        // if _linearize is true, S^2 terms are erased
      void add_term( const TERM_TYPE &_term, bool _linearize = false)
      {
        // exception: first term to be added to monome
        if ( terms.empty() )
        {
          terms.push_back( _term );
          return;
        }
      
        // positions new term -- terms is a sorted list
        typename std::list<TERM_TYPE> :: iterator i_where;
        i_where = std::find_if( terms.begin(), terms.end(), 
                                bind1st( std::less_equal<TERM_TYPE>(), _term ) );
       
        // insert new term in list
        if( i_where == terms.end() )
         terms.push_back(_term);
        else if ( _linearize and  *i_where == _term ) // linearizes here
          terms.erase(i_where);
        else
          terms.insert(i_where, _term);
      }

      // logic and algebraic operations
        // this function returns true if order of "this" is smaller than _comp
        // or if the variables in "this" are indexed smaller than those of _comp
      bool operator<(const Monome<_COEF, _TERM> & _comp) const
      {
        size_t mem_size = order();
        size_t _mem_size = _comp.order();
        
        if ( mem_size != _mem_size ) 
          return (mem_size < _mem_size);
      
        if ( mem_size == 0 ) // neither this nor _comp contain any terms
          return false;
      
        return std::lexicographical_compare( terms.begin(), terms.end(),
                                             _comp.terms.begin(),
                                             _comp.terms.end());
      }

      bool operator == (const Monome<_COEF, _TERM> &_comp) const
      {
        types::t_unsigned mem_size = order();
        types::t_unsigned _mem_size = _comp.order();
        if ( mem_size != _mem_size ) 
          return false;
        if ( mem_size == 0 )
          return true;
      
        return std::equal( terms.begin(), terms.end(),
                           _comp.terms.begin() );
      }
      bool operator<=(const Monome<_COEF, _TERM> & _comp) const
        { return ( operator< (_comp) || operator==(_comp) ); }
      void operator += (const Monome<_COEF, _TERM> &_comp)
        { coefficient += _comp.coefficient; }
      void operator = (const Monome<_COEF, _TERM> &_comp)
        { coefficient = _comp.coefficient; terms = _comp.terms; }
      void operator *= (const COEF_TYPE &factor)
        { coefficient *= factor; };
      void operator *= (const Monome<_COEF, _TERM> &factor)
        { coefficient *= factor; };
      void multiply(const Monome<_COEF, _TERM> &_mf, bool _linearize=false)
      {
        typename std::list<_TERM>:: const_iterator i_term = _mf.terms.begin();
        typename std::list<_TERM>:: const_iterator i_end = _mf.terms.end();
        for( ; i_term != i_end; ++i_term )
          add_term(*i_term, _linearize);
        coefficient *= _mf.coefficient;
      }

      // others
        // checks if monome contains _term
        // returns the number of instances _terms in monome
      bool contains_term (TERM_TYPE &_term)
      {
        if ( terms.empty() ) 
          return false;
      
        typename std::list<TERM_TYPE> :: iterator i_which;
        i_which = std::find( terms.begin(), terms.end(), _term);
      
        return (i_which != terms.end() );
      }
      types::t_int max_term() const
      {
       return  terms.empty() ? -1 : *(std::max_element(terms.begin(), terms.end(),
                                                       std::less<TERM_TYPE>()) );
      } 
      void print_terms( std::ostream &stream ) const
      {
        typename std::list <TERM_TYPE> :: const_iterator i_term = terms.begin();
        typename std::list <TERM_TYPE> :: const_iterator i_term_last = terms.end();
      
        for(; i_term != i_term_last; ++i_term )
          stream << " " << *i_term;
      }
      void print_out( std::ostream &stream ) const
      {
        if ( terms.empty() )
        {
          stream << " " << coefficient;
          return; 
        }
        typename std::list <TERM_TYPE> :: const_iterator i_term = terms.begin();
        typename std::list <TERM_TYPE> :: const_iterator i_term_last = terms.end();
      
        stream << " " << coefficient << "*";
        for(; i_term != i_term_last; i_term++ )
          stream << "S_" << *i_term;
      }
      std::ostream& operator<<( std::ostream &stream ) const
        { print_out(stream); return stream; }

  };
  

} // namespace function
#endif
