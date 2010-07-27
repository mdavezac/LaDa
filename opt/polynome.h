#ifndef _OPT_POLYNOME_H_
#define _OPT_POLYNOME_H_

#include "LaDaConfig.h"

#include <list>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cmath>
#include <math.h>


#include "function_base.h"
#include "monome.h"
#include <iomanip>

#include <opt/mpi.h>
#include <math/fuzzy.h>
#ifdef LADA_MPI
# include <opt/mpi.h>
# include <boost/lambda/lambda.hpp>
# include <boost/mpi/collectives.hpp>
# include <functional>
#endif 

namespace LaDa
{
  namespace function {

    //! \brief Defines a collection of monomials Monome and a polynomial function.
    //! \details This class provides two types of behaviors: 
    //!          - the behaviors expected from a collection of monomials as algebraic objects.
    //!          - the behaviors of a functional as defined by function::Base
    //!          .
    //!          When evaluating a polynomial, it is expected that the terms of
    //!          the monomials are indices to the variables from the
    //!          function::Base base class.
    //! \todo Separate the algebraic polynomial form the functional polynomial.
    //!       The second would be a derived class of the first. Requiring that
    //!       \a T_TERM is an index makes the templatization of T_TERM a bit
    //!       artificial.
    template<class T_TYPE=types::t_real, class T_TERM=types::t_int >
    class Polynome : public Base<T_TYPE> LADA_MPI_CODE( LADA_COMMA public MPI_COMMDEC )
    {
      protected:
        //! The type of the base class
        typedef Base<T_TYPE> t_Base;
        //! The type of this class
        typedef Polynome<T_TYPE, T_TERM> t_This;

      public:
        //! Type of the variables
        typedef T_TYPE t_Type;
        //! Type of each term
        typedef T_TERM t_Term;
        //! see function::Base::t_Container
        typedef typename Base<t_Type> :: t_Container t_Container;
        //! The type of the monomials
        typedef Monome<t_Type, t_Term> t_Monome;
        //! The type of the container of monomials
        typedef typename std::vector< t_Monome > t_Monomes;

      protected:
        using t_Base :: variables;

      protected: 
        //! The list of monomials
        t_Monomes monomes;

      public: 
        //! Constructor.
        Polynome() : t_Base() {}
        //! Constructor and variables initializer.
        Polynome   ( types::t_int nb )
                 : t_Base(nb) {}
        //! Copy constructor
        Polynome  ( const t_This &_p ) 
                 : monomes(_p.monomes)
                   LADA_MPI_CODE( LADA_COMMA MPI_COMMCOPY( _p ) ) {}

        //! Destructor
        virtual ~Polynome() {};
        
        //! Empties the container of monomials
        void clear();
        //! Returns a constant iterator to the first monomial.
        typename t_Monomes :: const_iterator monomes_begin() const { return monomes.begin(); } 
        //! Returns a constant iterator to the one past the last monomial.
        typename t_Monomes :: const_iterator monomes_end() const { return monomes.end(); } 

        //! Adds a monomial to the polynomial
        void  add(const t_Monome &_monome);
        //! Adds a collection of monomials to the polynomial.
        void  add(t_This &_polynome);
        //! Adds a collection of monomials to the polynomial with a minus sign.
        void  sub(t_This &_polynome);
        //! Adds a monomial to the polynomial
        void  operator += (const  t_Monome  &_monome) { add(_monome); } 
        //! Adds a collection of monomials to the polynomial.
        void  operator += (t_This &polynome) { add(polynome); }
        //! Adds a collection of monomials to the polynomial with a minus sign.
        void  operator -= (t_This &polynome)  { sub(polynome); }
        //! Multiplies all monomial coefficients by a factor.
        void  operator *= (const t_Type &factor);
        //! \brief Multiplies each monomial by the monomial \a _mf
        //! \details If \a _linearize is set to true, then the monomials are
        //!          linearized according to the description in the overview \a
        //!          opt::Monome.
        void  multiply(const t_Monome &_mf, bool _linearize=false);
        //! \brief Multiplies each monomial by the monomials of the polynomial \a _pf
        //! \details If \a _linearize is set to true, then the monomials are
        //!          linearized according to the description in the overview \a
        //!          opt::Monome.
        void  multiply(const t_This &_pf, bool _linearize=false);

        //! Evaluates the polynomial.
        virtual t_Type evaluate();
        //! Evaluates the gradient of the polynomial.
        virtual void evaluate_gradient(t_Type* const _i_grad);
        //! Evaluates the polynomial and its gradient.
        virtual t_Type evaluate_with_gradient( t_Type* const _i_grad);
        //! Evaluates the gradient of the polynomial in direction \a _pos
        virtual t_Type evaluate_one_gradient( types::t_unsigned _pos);

        //! Dumps the polynomial to a stream.
        void print_out( std::ostream &stream ) const;
        //! Dumps the polynomial to a stream.
        std::ostream& operator<<( std::ostream &stream ) const
          { print_out(stream); return stream; }
    };
   
    template<class T_TYPE, class T_TERM >
      inline void Polynome<T_TYPE, T_TERM> :: clear()
      {
        monomes.clear();
        if( variables ) variables->clear(); 
      }
    
    template<class T_TYPE, class T_TERM >
      inline void Polynome<T_TYPE, T_TERM> :: add(const t_Monome &_monome)
      {
        if ( monomes.empty() )
        {
          monomes.push_back(_monome);
          return;
        }
        
        typename t_Monomes :: iterator i_monome;
        i_monome = std::find_if( monomes.begin(), monomes.end(),
                                 bind1st( std::less_equal< t_Monome >(), _monome ) );
      
        if ( i_monome == monomes.end() ) // places _monome at end of sorted list
          monomes.push_back(_monome);           
        else if ( *i_monome == _monome )    // equivalent monome already exists
        {
          *i_monome += _monome;
          if ( math::eq(i_monome->coefficient, 0e0) )
            monomes.erase(i_monome);
        }
        else       // places _monome in the sorted list, before i_monome
          monomes.insert(i_monome,_monome);
      }
      
    template<class T_TYPE, class T_TERM >
      inline void Polynome<T_TYPE, T_TERM> :: add(t_This &_polynome)
      {
        typename t_Monomes :: iterator i_monome = _polynome.monomes.begin();
        typename t_Monomes :: iterator i_end = _polynome.monomes.end();
        for( ; i_monome != i_end; ++i_monome)
          add(*i_monome);
      }
    
    template<class T_TYPE, class T_TERM >
      inline void Polynome<T_TYPE, T_TERM> :: sub(t_This &_polynome)
      {
        typename t_Monomes :: iterator i_monome = _polynome.monomes.begin();
        typename std::list< t_Monome > :: iterator i_end = _polynome.monomes.end();
        for( ; i_monome != i_end; ++i_monome)
          { i_monome->coefficient *= -1.0; add(*i_monome); i_monome->coefficient *= -1.0; }
      }
      
    template<class T_TYPE, class T_TERM >
      inline void Polynome<T_TYPE, T_TERM> :: operator *= (const t_Type &factor)
      {
        typename t_Monomes :: iterator i_monome = monomes.begin();
        typename t_Monomes :: iterator i_end = monomes.end();
        for( ; i_monome != i_end; ++i_monome)
          i_monome->coefficient *= factor;
      }
      
    template<class T_TYPE, class T_TERM >
      inline void Polynome<T_TYPE, T_TERM> :: multiply(const t_Monome &_mf, bool _linearize)
      {
        t_This this_copy( *this );
        typename t_Monomes :: const_iterator i_monome = this_copy.monomes.begin();
        typename t_Monomes :: const_iterator i_end = this_copy.monomes.end();
        monomes.clear();
        for( ; i_monome != i_end; ++i_monome )
        {
          t_Monome monome( *i_monome );
          monome.multiply(_mf, _linearize);
          add(monome);
        }
      }
    
    template<class T_TYPE, class T_TERM >
      inline void Polynome<T_TYPE, T_TERM> :: multiply(const t_This &_pf, bool _linearize)
      {
        t_This this_copy( *this );
        typename t_Monomes :: const_iterator i_monome = _pf.monomes.begin();
        typename t_Monomes :: const_iterator i_end = _pf.monomes.end();
        typename t_Monomes :: const_iterator i_this_begin = this_copy.monomes.begin();
        typename t_Monomes :: const_iterator i_this_end = this_copy.monomes.end();
        typename t_Monomes :: const_iterator i_this;
        monomes.clear();
        for( ; i_monome != i_end; ++i_monome )
          for( i_this = i_this_begin; i_this != i_this_end; ++i_this )
          {
            t_Monome monome( *i_this );
            monome.multiply(*i_monome, _linearize);
            add(monome);
          }
      }

    template<class T_TYPE, class T_TERM >
      inline typename Polynome<T_TYPE, T_TERM> :: t_Type 
        Polynome<T_TYPE, T_TERM> :: evaluate() 
        {
          if ( monomes.empty() or not variables )
            return t_Type(0); 

          t_Type value(0);
        
          typename t_Monomes :: const_iterator i_monome = monomes.begin();
#         ifdef LADA_MPI
            typename t_Monomes :: const_iterator i_monome_end;
#         else 
            typename t_Monomes :: const_iterator i_monome_end = monomes.end();
#         endif
          typename t_Container :: const_iterator i_real = variables->begin();
          LADA_MPI_CODE(
            types :: t_unsigned nperproc = monomes.size() / MPI_COMM.size(); 
            types :: t_unsigned remainder = monomes.size() % MPI_COMM.size();
            i_monome +=  MPI_COMM.rank() * nperproc + std::min( (types::t_int) remainder,
                                                             MPI_COMM.rank() );
            i_monome_end = i_monome + nperproc;
            if( remainder and MPI_COMM.rank() < remainder ) ++i_monome_end;
          )
          for ( ; i_monome != i_monome_end; ++i_monome )
          {
            typename t_Monome :: t_Terms :: const_iterator i_term = i_monome->terms.begin();
            typename t_Monome :: t_Terms :: const_iterator i_term_end = i_monome->terms.end();
            t_Type v_monome = i_monome->coefficient;
        
            for( ; i_term != i_term_end; i_term++)
              v_monome *= t_Type( *(i_real + *i_term) );
        
            value += v_monome;
          } // end of loop over monomes
        
          LADA_MPI_CODE( value = boost::mpi::all_reduce( MPI_COMM, value,
                                                     std::plus<types::t_real>() ); )

          return value;
        }
      
    template<class T_TYPE, class T_TERM >
      inline void Polynome<T_TYPE, T_TERM> :: evaluate_gradient(t_Type* const _i_grad)
      {
        // clears gradient
        t_Type *ptr_grad = _i_grad;
        t_Type *ptr_end = _i_grad + variables->size();
        for(; ptr_grad != ptr_end; ptr_grad++)
          *ptr_grad = t_Type(0);

        if ( monomes.empty() or not variables ) return;
      
        t_Type g_monome;
         
        // loops over each monome.
#       ifdef LADA_ADDHERE
#         error Please change LADA_ADDHERE to something not already defined
#       endif
        typename t_Monomes :: const_iterator i_monome = monomes.begin();
#       ifdef LADA_MPI
          typename t_Monomes :: const_iterator i_monome_end;
#         define LADA_ADDHERE array
#       else 
          typename t_Monomes :: const_iterator i_monome_end = monomes.end();
#         define LADA_ADDHERE _i_grad
#       endif
        typename t_Container :: const_iterator i_real = variables->begin();
        LADA_MPI_CODE(
          types::t_unsigned nperproc = monomes.size() / MPI_COMM.size(); 
          types :: t_unsigned remainder = monomes.size() % MPI_COMM.size();
          i_monome +=  MPI_COMM.rank() * nperproc + std::min( (types::t_int) remainder,
                                                           MPI_COMM.rank() );
          i_monome_end = i_monome + nperproc;
          if( remainder and MPI_COMM.rank() < remainder ) ++i_monome_end;
          types::t_real *LADA_ADDHERE = new types::t_real[ monomes.size() ];
          std::fill( LADA_ADDHERE, LADA_ADDHERE + monomes.size(), types::t_real(0) );
        )
        for ( ; i_monome != i_monome_end; ++i_monome ) // monome loop 
        {
          typename t_Monome :: t_Terms :: const_iterator i_term_begin = i_monome->terms.begin();
          typename t_Monome :: t_Terms :: const_iterator i_term_end = i_monome->terms.end();
          typename t_Monome :: t_Terms :: const_iterator i_excluded = i_term_begin;
          typename t_Monome :: t_Terms :: const_iterator i_included;
      
        
          // outer loop runs over derivatives
          for( i_excluded = i_term_begin; i_excluded != i_term_end;
               ++i_excluded)
          {
            // sets monome gradient storage
            g_monome = i_monome->coefficient;
      
            // inner loop runs over other variables in monome
            for( i_included = i_term_begin; i_included != i_term_end;
                 ++i_included )
              if ( i_included != i_excluded )
                g_monome *= t_Type( *(i_real + *i_included) );
      
            // now sums gradient
            *(LADA_ADDHERE + *i_excluded) += g_monome;
      
          } // end of outer "derivative" loop
        } // end of loop over monomes
        LADA_MPI_CODE
        ( 
          // boost does not do in-place reduction.
          BOOST_MPI_CHECK_RESULT
          (
            MPI_Allreduce,
            (
              MPI_IN_PLACE, LADA_ADDHERE, monomes.size(),
              boost::mpi::get_mpi_datatype<types::t_real>(*LADA_ADDHERE),
              boost::mpi::is_mpi_op< std::plus<types::t_real>, types::t_real>::op(),
              (MPI_Comm) MPI_COMM
            )
          );
          std::transform( _i_grad, _i_grad + monomes.size(), LADA_ADDHERE, _i_grad,
                          boost::lambda::_1 + boost::lambda::_2 ); 
          delete[] LADA_ADDHERE;
        )
#       undef LADA_ADDHERE
      }
    
    template<class T_TYPE, class T_TERM >
      inline typename Polynome<T_TYPE, T_TERM> :: t_Type 
        Polynome<T_TYPE, T_TERM> :: evaluate_with_gradient( t_Type* const _i_grad) 
        {
          // clears gradient
          t_Type *ptr_grad = _i_grad;
          t_Type *ptr_end = _i_grad + variables->size();
          for(; ptr_grad != ptr_end; ptr_grad++)
            *ptr_grad = t_Type(0);

          if ( monomes.empty() or not variables )
            return t_Type(0);
        
          t_Type v_monome, g_monome;
          t_Type value(0);
           
          // loops over each monome.
          typename t_Monomes :: const_iterator i_monome = monomes.begin();
#         ifdef LADA_MPI
            typename t_Monomes :: const_iterator i_monome_end;
#           define LADA_ADDHERE array
#         else 
            typename t_Monomes :: const_iterator i_monome_end = monomes.end();
#           define LADA_ADDHERE _i_grad
#         endif
          typename t_Container :: const_iterator i_real = variables->begin();
          LADA_MPI_CODE(
            types :: t_unsigned nperproc = monomes.size() / MPI_COMM.size(); 
            types :: t_unsigned remainder = monomes.size() % MPI_COMM.size();
            i_monome +=  MPI_COMM.rank() * nperproc + std::min( (types::t_int) remainder,
                                                             MPI_COMM.rank() );
            i_monome_end = i_monome + nperproc;
            if( remainder and MPI_COMM.rank() < remainder ) ++i_monome_end;
            types::t_real *LADA_ADDHERE = new types::t_real[ monomes.size() ];
          )
          for ( ; i_monome != i_monome_end; i_monome++ ) // monome loop 
          {
            typename t_Monome :: t_Terms :: const_iterator i_term_begin = i_monome->terms.begin();
            typename t_Monome :: t_Terms :: const_iterator i_term_end = i_monome->terms.end();
            typename t_Monome :: t_Terms :: const_iterator i_excluded = i_term_begin;
            typename t_Monome :: t_Terms :: const_iterator i_included;
        
            v_monome = i_monome->coefficient;
            
            // outer loop runs over derivatives
            for(i_excluded = i_term_begin; i_excluded!=i_term_end; ++i_excluded)
            {
              // sets monome gradient storage
              g_monome = i_monome->coefficient;
        
              // inner loop runs over other variables in monome
              for( i_included = i_term_begin; i_included!=i_term_end; ++i_included )
                if ( i_included != i_excluded )
                  g_monome *= t_Type( *(i_real + *i_included) );
        
              // now sums gradient
              *(_i_grad + *i_excluded) += g_monome;
        
              // evaluate monome
              v_monome *= t_Type( *(i_real + *i_excluded) );
            }
        
            value += v_monome;
        
          } // end of loop over monomes
        
          LADA_MPI_CODE
          ( 
            BOOST_MPI_CHECK_RESULT
            (
              MPI_Allreduce,
              (
                MPI_IN_PLACE, LADA_ADDHERE, monomes.size(),
                boost::mpi::get_mpi_datatype<types::t_real>(*LADA_ADDHERE),
                boost::mpi::is_mpi_op< std::plus<types::t_real>, types::t_real>::op(),
                (MPI_Comm) MPI_COMM
              )
            );
            std::transform( _i_grad, _i_grad + monomes.size(), LADA_ADDHERE, _i_grad,
                            boost::lambda::_1 + boost::lambda::_2 ); 
            delete[] LADA_ADDHERE;
            value = boost::mpi::all_reduce( MPI_COMM, value,
                                            std::plus<types::t_real>() ); 
          )
#         undef LADA_ADDHERE
          return value;
        }
        
    template<class T_TYPE, class T_TERM >
      inline typename Polynome<T_TYPE, T_TERM> :: t_Type 
        Polynome<T_TYPE, T_TERM> :: evaluate_one_gradient( types::t_unsigned _pos)
        {
          // clears gradient
          t_Type result = t_Type(0);

          if ( monomes.empty() or not variables ) return t_Type(0);
        
           
          // loops over each monome.
          typename t_Monomes :: const_iterator i_monome = monomes.begin();
#         ifdef LADA_MPI
            typename t_Monomes :: const_iterator i_monome_end;
#         else 
            typename t_Monomes :: const_iterator i_monome_end = monomes.end();
#         endif
          typename t_Container :: const_iterator i_real = variables->begin();
          LADA_MPI_CODE(
            types :: t_unsigned nperproc = monomes.size() / MPI_COMM.size(); 
            types :: t_unsigned remainder = monomes.size() % MPI_COMM.size();
            i_monome +=  MPI_COMM.rank() * nperproc + std::min( (types::t_int) remainder,
                                                             MPI_COMM.rank() );
            i_monome_end = i_monome + nperproc;
            if( remainder and MPI_COMM.rank() < remainder ) ++i_monome_end;
          )
          for ( ; i_monome != i_monome_end; i_monome++ ) // monome loop 
          {
            typename t_Monome :: t_Terms :: const_iterator i_term = i_monome->terms.begin();
            typename t_Monome :: t_Terms :: const_iterator i_term_end = i_monome->terms.end();
            t_Type partial_grad = t_Type(0);
        
            t_Term exponent = std::count( i_term, i_term_end, t_Term(_pos) );
            
            if ( exponent > 0 )
            {
              partial_grad = i_monome->coefficient;
              i_term = i_monome->terms.begin();
              for( ; i_term != i_term_end; ++i_term )
                if ( *i_term != t_Term(_pos) )
                  partial_grad *= t_Type( *(i_real + *i_term ) );

              partial_grad *= t_Type(exponent);
              t_Type value = *(i_real + _pos);
              for( t_Term n=1; n < exponent; ++i_term )
                partial_grad *= value;
            }
        
            result += partial_grad;
        
          } // end of loop over monomes
          LADA_MPI_CODE( result = boost::mpi::all_reduce( MPI_COMM, result,
                                                      std::plus<types::t_real>() ); )
        
          return result;
        }

    template<class T_TYPE, class T_TERM >
      inline void Polynome<T_TYPE, T_TERM> :: print_out( std::ostream &stream ) const
      {
        typename t_Monomes :: const_iterator i_monome = monomes.begin();
        typename t_Monomes :: const_iterator i_monome_end = monomes.end();
        types::t_int i = 1;
      
        std::cout << std::right << std::noshowpos;
        stream << "Polynome: ";
        i_monome->print_out(stream);
        for(++i_monome; i_monome != i_monome_end; ++i_monome )
        {
          stream << " +";
          i_monome->print_out(stream);
          i++;
          if ( i % 3 == 0 )
            stream << std::endl << "            ";
        }  
        stream << std::endl;
      }
      
  } // namespace opt
} // namespace Lada
#endif
