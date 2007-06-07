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
#ifndef _OPT_POLYNOME_H_
#define _OPT_POLYNOME_H_

#include <list>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cmath>
#include <math.h>
#include <opt/opt_function_base.h>
#include <opt/opt_monome.h>

#include <mpi/mpi_object.h>

namespace function {

  template<class _TYPE=types::t_real, class _TERM=types::t_int >
  class Polynome : public Base<_TYPE>
  {
    public:
      typedef _TYPE t_Type;
      typedef _TERM t_Term;
      typedef typename Base<t_Type> :: t_Container t_Container;
      typedef typename std::list<Monome<t_Type, t_Term> > Monome_Container;
      typedef typename std::list<Monome<t_Type, t_Term> > :: iterator Monome_Iterator;
      typedef typename std::list<Monome<t_Type, t_Term> > :: const_iterator const_Monome_Iterator;

    protected:
      using Base<t_Type> :: variables;

    protected: 
      Monome_Container monomes;

    public: 
      Polynome() : Base<t_Type>() {};
      Polynome( types::t_int nb ) : Base<t_Type>(nb) {};
      Polynome( const Polynome<t_Type, t_Term> &_p ) : monomes(_p.monomes) {};
      virtual ~Polynome() {};
      
      void operator=(const Polynome<t_Type, t_Term> &_p ) 
      {
        monomes.clear();
        monomes = _p.monomes;
      }


      // setting things up
      void clear()
      {
        monomes.clear();
        if( variables ) 
          variables->clear(); 
      }
      const_Monome_Iterator monomes_begin() const
        { return monomes.begin(); } 
      const_Monome_Iterator monomes_end() const
        { return monomes.end(); } 

      // algebraic operations
      void  add(const Monome<t_Type, t_Term> &_add)
      {
        if ( monomes.empty() )
        {
          monomes.push_back(_add);
          return;
        }
        
        typename std::list< Monome<t_Type, t_Term> > :: iterator i_monome;
        i_monome = std::find_if( monomes.begin(), monomes.end(),
                                 bind1st( std::less_equal< Monome<t_Type, t_Term> >(), _add ) );
      
        if ( i_monome == monomes.end() ) // places _add at end of sorted list
          monomes.push_back(_add);           
        else if ( *i_monome == _add )    // equivalent monome already exists
        {
          *i_monome += _add;
          if ( is_null(i_monome->coefficient) )
            monomes.erase(i_monome);
        }
        else       // places _add in the sorted list, before i_monome
          monomes.insert(i_monome,_add);
      }
      void  add(Polynome<t_Type, t_Term> &_add)
      {
        typename std::list< Monome<t_Type, t_Term> > :: iterator i_monome = _add.monomes.begin();
        typename std::list< Monome<t_Type, t_Term> > :: iterator i_last = _add.monomes.end();
        for( ; i_monome != i_last; ++i_monome)
          add(*i_monome);
      }
      void  sub(Polynome<t_Type, t_Term> &_add)
      {
        typename std::list< Monome<t_Type, t_Term> > :: iterator i_monome = _add.monomes.begin();
        typename std::list< Monome<t_Type, t_Term> > :: iterator i_last = _add.monomes.end();
        for( ; i_monome != i_last; ++i_monome)
          { i_monome->coefficient *= -1.0; add(*i_monome); i_monome->coefficient *= -1.0; }
      }
      void  operator += (const  Monome<t_Type, t_Term>  &_monome)
        { add(_monome); } 
      void  operator += (Polynome<t_Type, t_Term> &polynome)
        { add(polynome); }
      void  operator -= (Polynome<t_Type, t_Term> &polynome)
        { sub(polynome); }
      void  operator *= (const t_Type &factor)
      {
        typename std::list<Monome<t_Type, t_Term> >:: iterator i_monome = monomes.begin();
        typename std::list<Monome<t_Type, t_Term> >:: iterator i_end = monomes.end();
        for( ; i_monome != i_end; ++i_monome)
          i_monome->coefficient *= factor;
      }
      void  multiply(const Monome<t_Type, t_Term> &_mf, bool _linearize=false)
      {
        Polynome<t_Type, t_Term> this_copy( *this );
        typename std::list<Monome<t_Type, t_Term> >:: const_iterator i_monome = this_copy.monomes.begin();
        typename std::list<Monome<t_Type, t_Term> >:: const_iterator i_end = this_copy.monomes.end();
        monomes.clear();
        for( ; i_monome != i_end; ++i_monome )
        {
          Monome< t_Type, t_Term> monome( *i_monome );
          monome.multiply(_mf, _linearize);
          add(monome);
        }
      }
      void  multiply(const Polynome<t_Type, t_Term> &_pf, bool _linearize=false)
      {
        Polynome<t_Type, t_Term> this_copy( *this );
        typename std::list<Monome<t_Type, t_Term> >:: const_iterator i_monome = _pf.monomes.begin();
        typename std::list<Monome<t_Type, t_Term> >:: const_iterator i_end = _pf.monomes.end();
        typename std::list<Monome<t_Type, t_Term> >:: const_iterator i_this_begin = this_copy.monomes.begin();
        typename std::list<Monome<t_Type, t_Term> >:: const_iterator i_this_end = this_copy.monomes.end();
        typename std::list<Monome<t_Type, t_Term> >:: const_iterator i_this;
        monomes.clear();
        for( ; i_monome != i_end; ++i_monome )
          for( i_this = i_this_begin; i_this != i_this_end; ++i_this )
          {
            Monome< t_Type, t_Term> monome( *i_this );
            monome.multiply(*i_monome, _linearize);
            add(monome);
          }
      }

      // evaluations
      virtual t_Type evaluate() // returns result
      {
        if ( monomes.empty() or not variables )
          return t_Type(0); 

        t_Type value(0);
      
        typename std::list< Monome<t_Type, t_Term> > :: const_iterator i_monome = monomes.begin();
        typename std::list< Monome<t_Type, t_Term> > :: const_iterator i_monome_last = monomes.end();
        typename t_Container :: const_iterator i_real = variables->begin();
        for ( ; i_monome != i_monome_last; i_monome++ )
        {
          typename std::list<t_Term> :: const_iterator i_term = i_monome->terms.begin();
          typename std::list<t_Term> :: const_iterator i_term_last = i_monome->terms.end();
          t_Type v_monome = i_monome->coefficient;
      
          for( ; i_term != i_term_last; i_term++)
            v_monome *= t_Type( *(i_real + *i_term) );
      
          value += v_monome;
        } // end of loop over monomes
      
        return value;
      }
      virtual void evaluate_gradient(t_Type* const _i_grad) // result returns in gradient
      {
        // clears gradient
        t_Type *ptr_grad = _i_grad;
        t_Type *ptr_last = _i_grad + variables->size();
        for(; ptr_grad != ptr_last; ptr_grad++)
          *ptr_grad = t_Type(0);

        if ( monomes.empty() or not variables )
          return;
      
        t_Type g_monome;
         
        // loops over each monome.
        typename std::list< Monome<t_Type, t_Term> > :: const_iterator i_monome = monomes.begin();
        typename std::list< Monome<t_Type, t_Term> > :: const_iterator i_monome_last = monomes.end();
        typename t_Container :: const_iterator i_real = variables->begin();
        for ( ; i_monome != i_monome_last; ++i_monome ) // monome loop 
        {
          typename std::list< t_Term> :: const_iterator i_term_begin = i_monome->terms.begin();
          typename std::list< t_Term> :: const_iterator i_term_last = i_monome->terms.end();
          typename std::list< t_Term> :: const_iterator i_excluded = i_term_begin;
          typename std::list< t_Term> :: const_iterator i_included;
      
        
          // outer loop runs over derivatives
          for( i_excluded = i_term_begin; i_excluded != i_term_last;
               ++i_excluded)
          {
            // sets monome gradient storage
            g_monome = i_monome->coefficient;
      
            // inner loop runs over other variables in monome
            for( i_included = i_term_begin; i_included != i_term_last;
                 ++i_included )
              if ( i_included != i_excluded )
                g_monome *= t_Type( *(i_real + *i_included) );
      
            // now sums gradient
            *(_i_grad + *i_excluded) += g_monome;
      
          } // end of outer "derivative" loop
        } // end of loop over monomes
      }
      // returns evaluation, and stores gradient directions in gradient
      virtual t_Type evaluate_with_gradient( t_Type* const _i_grad) 
      {
        // clears gradient
        t_Type *ptr_grad = _i_grad;
        t_Type *ptr_last = _i_grad + variables->size();
        for(; ptr_grad != ptr_last; ptr_grad++)
          *ptr_grad = t_Type(0);

        if ( monomes.empty() or not variables )
          return t_Type(0);
      
        t_Type v_monome, g_monome;
        t_Type value(0);
         
        // loops over each monome.
        typename std::list< Monome<t_Type, t_Term> > :: const_iterator i_monome = monomes.begin();
        typename std::list< Monome<t_Type, t_Term> > :: const_iterator i_monome_last = monomes.end();
        typename t_Container :: const_iterator i_real = variables->begin();
        for ( ; i_monome != i_monome_last; i_monome++ ) // monome loop 
        {
          typename std::list<t_Term> :: const_iterator i_term_begin = i_monome->terms.begin();
          typename std::list<t_Term> :: const_iterator i_term_last = i_monome->terms.end();
          typename std::list<t_Term> :: const_iterator i_excluded = i_term_begin;
          typename std::list<t_Term> :: const_iterator i_included;
      
          v_monome = i_monome->coefficient;
          
          // outer loop runs over derivatives
          for(i_excluded = i_term_begin; i_excluded!=i_term_last; ++i_excluded)
          {
            // sets monome gradient storage
            g_monome = i_monome->coefficient;
      
            // inner loop runs over other variables in monome
            for( i_included = i_term_begin; i_included!=i_term_last; ++i_included )
              if ( i_included != i_excluded )
                g_monome *= t_Type( *(i_real + *i_included) );
      
            // now sums gradient
            *(_i_grad + *i_excluded) += g_monome;
      
            // evaluate monome
            v_monome *= t_Type( *(i_real + *i_excluded) );
          }
      
          value += v_monome;
      
        } // end of loop over monomes
      
        return value;
      }
      
      virtual t_Type evaluate_one_gradient( types::t_unsigned _pos)
      {
        // clears gradient
        t_Type result;

        if ( monomes.empty() or not variables )
          return t_Type(0);
      
         
        // loops over each monome.
        typename std::list< Monome<t_Type, t_Term> > :: const_iterator i_monome = monomes.begin();
        typename std::list< Monome<t_Type, t_Term> > :: const_iterator i_monome_last = monomes.end();
        typename t_Container :: const_iterator i_real = variables->begin();
        for ( ; i_monome != i_monome_last; i_monome++ ) // monome loop 
        {
          typename std::list<t_Term> :: const_iterator i_term = i_monome->terms.begin();
          typename std::list<t_Term> :: const_iterator i_term_last = i_monome->terms.end();
          t_Type partial_grad;
      
          t_Term exponent = std::count( i_term, i_term_last, t_Term(_pos) );
          
          if ( exponent > 0 )
          {
            partial_grad = i_monome->coefficient;
            i_term = i_monome->terms.begin();
            for( ; i_term != i_term_last; ++i_term )
              if ( *i_term != t_Term(_pos) )
                partial_grad *= t_Type( *(i_real + *i_term ) );

            partial_grad *= t_Type(exponent);
            t_Type value = *(i_real + _pos);
            for( t_Term n=1; n < exponent; ++i_term )
              partial_grad *= value;
          }
      
          result += partial_grad;
      
        } // end of loop over monomes
      
        return result;
      }

      // other
      void print_out( std::ostream &stream ) const
      {
        typename std::list < Monome<t_Type, t_Term> > :: const_iterator i_monome = monomes.begin();
        typename std::list < Monome<t_Type, t_Term> > :: const_iterator i_monome_last = monomes.end();
        types::t_int i = 1;
      
        std::cout << std::right << std::noshowpos;
        stream << "Polynome: ";
        i_monome->print_out(stream);
        for(++i_monome; i_monome != i_monome_last; ++i_monome )
        {
          stream << " +";
          i_monome->print_out(stream);
          i++;
          if ( i % 3 == 0 )
            stream << std::endl << "            ";
        }  
//       stream << " = " << evaluate();
        stream << std::endl;
      }
      std::ostream& operator<<( std::ostream &stream ) const
        { print_out(stream); return stream; }

    protected:
      bool is_null (t_Type& _coef) const
      {
        return std::abs(_coef) < 1e-6;
      }

    public: // MPI stuff
     bool broadcast( mpi::BroadCast &_bc )
     {
#ifdef _MPI
       types::t_int n = monomes.size();
       if ( not _bc.serialize( n ) ) return false;
       if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
         monomes.resize(n);
       Monome_Iterator i_monome = monomes.begin();
       Monome_Iterator i_monome_end = monomes.end();
       for(; i_monome != i_monome_end; ++i_monome )
         if(    ( not _bc.serialize( i_monome->coefficient ) )
             or ( not _bc.serialize( i_monome->terms ) ) ) return false;
#endif 
       return true;
     }
  };
 
} // namespace opt
#endif
