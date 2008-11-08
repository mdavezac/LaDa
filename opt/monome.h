//
//  Version: $Id$
//
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

namespace LaDa
{
  namespace function {

    template<class T_TYPE, class T_TERM> class Polynome;

    /** \brief Defines a Monomial.
     *  \details A monomial consists a coefficient and an array of terms. In
     *           general, the terms are indices meant to refer to the variables of a
     *           polynomial or higher order %function incorporating Monomes. They
     *           can be most anything though with default ordering and comparison
     *           operator, as well as a copy constructor.
     *   
     *           The monomials have the ability to linearize themselves, with the
     *           assumption that for any term \e a, we have \f$a^{2n+1}=a\f$ and
     *           \f$a^{2n} = 1\f$. Say you have a monomial containaing terms \e A,
     *           \e B, \e C, and \e A again. Then the linearized monomial would
     *           containe only \e B and \e C. Basically, think of a collection of
     *           +/-1 spins... */
    template<class T_COEF=types::t_real, class T_TERM=types::t_int >
    class Monome
    {
      friend class Polynome<T_COEF, T_TERM>;

      public:
        //! Type of the coefficient
        typedef T_COEF t_Coef;
        //! Type if the terms
        typedef T_TERM t_Term;
        //! Type of  a collection of indices
        typedef typename std::list<T_TERM> t_Terms;

      protected:
        //! Type of this class.
        typedef Monome<t_Coef, t_Term> t_This;

      protected:
        //! The factor of  the monomial
        t_Coef coefficient;
        //! The list of terms
        std::list<t_Term> terms;

      public:
        //! Constructor
        Monome() : coefficient(0) {};
        //! Constructor and Initializer
        Monome(const t_Coef &_coef) : coefficient(_coef) {};
        //! Copy Constructor.
        Monome(const t_This &_m) : coefficient (_m.coefficient), terms(_m.terms) {};
        //! Destructor
        ~Monome(){};

        //! Sets the coefficient to \a _c
        void set_coefficient( t_Coef _c )  { coefficient = _c; };
        //! Returns the value of the coefficient
        t_Coef get_coefficient() const { return coefficient; };
        //! Returns the order of the monomial
        size_t order() const { return terms.empty() ? 0 : terms.size(); }
        //! Returns a constant iterator to the first term.
        typename t_Terms :: const_iterator terms_begin() const { return terms.begin(); }
        //! Returns an iterator to the first term.
        typename t_Terms :: iterator terms_begin()   { return terms.begin(); }
        //! Returns a constant iterator to one past the last term.
        typename t_Terms :: const_iterator terms_end() const { return terms.end(); }
        //! Returns an iterator to one past the last term.
        typename t_Terms :: iterator terms_end() { return terms.end(); }

        //! \brief Adds a term to to the monome. 
        //! \details this is an index to the variables defined in a Polynome object
        //!          first add_term checks for memory availability
        //!          then it places the value of _term in a sorted list
        //!          if _linearize is true, S^2 terms are erased as explained in
        //!          opt::Monome overview.
        void add_term( const t_Term &_term, bool _linearize = false);

        //! \brief Strict weak ordering
        //! \details This function returns true if the order of this monomial is
        //!          strictly smaller than the order of \a _comp. If the orders
        //!          are equal, then a lexicographic comparison of the terms is
        //!          returned. 
        bool operator<(const t_This & _comp) const;

        //! \brief Order comparison followed by term comparison.
        //! \details Two monomials are equal if and only if they are of equal
        //!          order and their terms are equal.
        bool operator == (const t_This &_comp) const;

        //! \brief Uses Monome::operator<() and Monome::operator==()
        bool operator<=(const t_This & _comp) const
          { return ( operator< (_comp) || operator==(_comp) ); }
        //! \brief Sums the coefficients of monomial \a _comp into the
        //!        coefficient of this monomial.
        //! \warning The coefficients of two monomials can be summed if and only
        //!          if the two monomials are equal. This function does the
        //!          coefficient summing regardless...
        void operator += (const t_This &_comp) { coefficient += _comp.coefficient; }
        //! \brief Multiplies the coefficient by \a _factor
        void operator *= (const t_Coef &_factor) { coefficient *= _factor; };
        //! \brief Multiplies the coefficient by the coefficient of \a _factor
        void operator *= (const t_This &_factor) { coefficient *= _factor; };
        //! \brief Multiplies the full monomial by \a _mf
        //! \details In practice, this multiplying one coefficient by the other
        //!          and adding the terms of \a _mf to those already here. As a
        //!          bonus, if \a _linearize is set to true, then the resulting
        //!          monomial is linearized, in the sense explained in opt::Monome.
        void multiply(const t_This &_mf, bool _linearize=false);

        //! Returns true if the monomial contains \a _term
        bool contains_term (t_Term &_term);

        //! \brief Returns the "maximum" term.
        //! \details Maximum depends on the default ordering operator< for the
        //!          term.
        types::t_int max_term() const;

        //! Dumps the terms of the monomial to a stream.
        void print_terms( std::ostream &stream ) const;
        //! Dumps the monomial to a stream.
        void print_out( std::ostream &stream ) const;
        //! Dumps the monomial to a stream.
        std::ostream& operator<<( std::ostream &stream ) const
          { print_out(stream); return stream; }

    };
    
    template<class T_COEF, class T_TERM >
      inline void Monome<T_COEF, T_TERM> :: add_term( const t_Term &_term,
                                                      bool _linearize )
      {
        // exception: first term to be added to monome
        if ( terms.empty() )
        {
          terms.push_back( _term );
          return;
        }
      
        // positions new term -- terms is a sorted list
        typename t_Terms :: iterator i_where;
        i_where = std::find_if( terms.begin(), terms.end(), 
                                bind1st( std::less_equal<t_Term>(), _term ) );
       
        // insert new term in list
        if( i_where == terms.end() )
         terms.push_back(_term);
        else if ( _linearize and  *i_where == _term ) // linearizes here
          terms.erase(i_where);
        else
          terms.insert(i_where, _term);
      }

    template<class T_COEF, class T_TERM >
      inline bool Monome<T_COEF, T_TERM> :: operator<(const t_This & _comp) const
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

    template<class T_COEF, class T_TERM >
      inline bool Monome<T_COEF, T_TERM> :: operator==(const t_This &_comp) const
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
        
    template<class T_COEF, class T_TERM >
      inline void Monome<T_COEF, T_TERM> :: multiply(const t_This &_mf, 
                                                     bool _linearize   )
      {
        typename t_Terms:: const_iterator i_term = _mf.terms.begin();
        typename t_Terms:: const_iterator i_end = _mf.terms.end();
        for( ; i_term != i_end; ++i_term )
          add_term(*i_term, _linearize);
        coefficient *= _mf.coefficient;
      }

        
    template<class T_COEF, class T_TERM >
      inline bool Monome<T_COEF, T_TERM> :: contains_term (t_Term &_term)
      {
        if ( terms.empty() ) 
          return false;
      
        typename t_Terms :: iterator i_which;
        i_which = std::find( terms.begin(), terms.end(), _term);
      
        return (i_which != terms.end() );
      }

    template<class T_COEF, class T_TERM >
      inline types::t_int Monome<T_COEF, T_TERM> :: max_term() const
      {
        if ( terms.empty() ) return -1;
        return *( std::max_element(terms.begin(), terms.end(),
                                   std::less<t_Term>())        );
      } 
    
    template<class T_COEF, class T_TERM >
      inline void Monome<T_COEF, T_TERM> :: print_terms( std::ostream &stream ) const
      {
        typename t_Terms :: const_iterator i_term = terms.begin();
        typename t_Terms :: const_iterator i_term_last = terms.end();
      
        for(; i_term != i_term_last; ++i_term )
          stream << " " << *i_term;
      }
      
    template<class T_COEF, class T_TERM >
      inline void Monome<T_COEF, T_TERM> :: print_out( std::ostream &stream ) const
      {
        if ( terms.empty() )
        {
          stream << " " << coefficient;
          return; 
        }
        typename t_Terms :: const_iterator i_term = terms.begin();
        typename t_Terms :: const_iterator i_term_last = terms.end();
      
        stream << " " << coefficient << "*";
        for(; i_term != i_term_last; i_term++ )
          stream << "S_" << *i_term;
      }

  } // namespace function
} // namespace LaDa
#endif
