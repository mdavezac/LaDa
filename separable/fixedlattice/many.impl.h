//
//  Version: $Id$
//
#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/matrix.hpp>
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/numeric/ublas/operation.hpp>
#include<boost/regex.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/tokenizer.hpp>
#include<boost/fusion/algorithm/iteration/for_each.hpp>
#include<boost/spirit/home/phoenix/core.hpp>
#include<boost/spirit/home/phoenix/stl/algorithm/iteration.hpp>
#include<boost/spirit/home/phoenix/scope.hpp>
#include<boost/spirit/home/phoenix/container.hpp>
#include<boost/spirit/home/phoenix/operator.hpp>
#include<boost/spirit/home/phoenix/statement/if.hpp>

#include "prepare.h"

namespace CE
{
    INMANY2(template< class T_MATRIX, class T_VECTOR > void) 
      :: create_A_n_b( T_MATRIX &_A, T_VECTOR &_b )
      {
        namespace bf = boost::fusion;
        namespace bblas = boost::numeric::ublas;
        __DEBUGTRYBEGIN
        // Loop over inequivalent configurations.
        t_Vector X( current_dof() );
        std::fill( _A.data().begin(), _A.data().end(), 
                   typename t_Matrix::value_type(0) );
        std::fill( _b.data().begin(), _b.data().end(), 
                   typename t_Vector::value_type(0) );
        for( size_t i(0); i < mapping().size(); ++i )
        {
          // allows leave-one-out, or leave-many-out.
          if( mapping().do_skip(i) )  continue;

          bf::transform( *collapses_, ApplyCreateAnB<t_Vector>( X, i, *this ) );
          
          _A += mapping().weight(i) * bblas::outer_prod( X, X ); 
          _b += mapping().weight(i) * mapping().target(i) * X;
        }
        __DEBUGTRYEND(, "Error in Many::create_A_n_b()\n" )
      }

    INMANY2( template< class T_MATRIX, class T_VECTOR > void )
      :: operator()( T_MATRIX &_A, T_VECTOR &_b, types::t_unsigned _dim )
      {
        __DEBUGTRYBEGIN
        namespace bblas = boost::numeric::ublas;
        dim = _dim;
//       boost::fusion::for_each( collapses_, U_MFUNC( assign )(dim) );
        create_A_n_b( _A, _b );

        ApplyRegularization< T_MATRIX, T_VECTOR > applyreg( _A, _b, *this );
        boost::fusion::transform( *collapses_, applyreg );
        __DEBUGTRYEND(, "Error in Many::operator()()\n" )
      }
    
    INMANY( opt::ErrorTuple ) :: evaluate() const 
    {
      __DEBUGTRYBEGIN
      opt::ErrorTuple error;
      for(size_t n(0); n < mapping().size(); ++n )
        if( not mapping().do_skip(n) )
          error += opt::ErrorTuple( mapping().target( n ) - evaluate(n),
                                    mapping().weight( n ) );
      return error;
      __DEBUGTRYEND(, "Error in Many::evaluate()\n" )
    }

    INMANY( template< size_t _index > size_t ) ::  addone()
    {
      namespace bm = boost::mpl;
      namespace bf = boost::fusion;
      typedef boost::mpl::int_< _index > N;
      typedef typename bm::at<t_Collapses, N> :: type t_collapse;
      typedef typename bm::at<t_Separables, N> :: type t_separables;

      bf::at<N>( *collapses_ ).push_back( new typename bm::at<t_Collapses, N> :: type );
      bf::at<N>( *separables_ ).push_back( new t_separables );
      bf::at<N>( *collapses_ ).back().init( bf::at<N>( *separables_ ).back() );
      return bf::at<N>( *separables_ ).size() - 1;
    }

    INMANY( template< size_t _index > void )
      :: init( size_t _rank, size_t _dimensions )
      {
        namespace bblas = boost::numeric::ublas;
        namespace bf = boost::fusion;
        namespace bm = boost::mpl;
        // last collapse is set to fake range fake range 
        const bblas :: range a( 0, _rank );
        const bblas :: range b( 0, _dimensions );
        typedef boost::mpl::int_< _index > N;
        typedef typename bm::at<t_Collapses, N> :: type t_collapse;
        t_collapse &collapse( bf::at<N>( *collapses_).back() );
        collapse.coefficients_interface().set( coefficients_, a, b );
      
        // Now computes rank and dimensions.
        coefficients_.resize( dof(), dimensions() );
        boost::fusion::transform( *collapses_,
                                  ApplyResize<const t_Coefficients>( coefficients_ ) );
      }
    INMANY( typename MANYHEAD::t_Matrix::value_type ) :: evaluate( size_t _n) const 
    {
      namespace bp = boost::phoenix;
      namespace bpa = boost::phoenix::arg_names;
      namespace bpl = boost::phoenix::local_names;
      bp::function< PHOENIX_MEMFUNC1( evaluate ) > memfunc_eval;
      types::t_real result(0);
      boost::fusion::for_each
      ( 
        const_viewasref(*collapses_),
        bp::for_each
        ( 
          bpa::arg1,
          bp::lambda[ bp::ref(result) += memfunc_eval( bpa::arg1, bp::ref(_n) ) ]
        ) 
      );
      return result;
    }
    INMANY( size_t ) :: current_dof() const 
    {
      namespace bp = boost::phoenix;
      namespace bpa = boost::phoenix::arg_names;
      namespace bpl = boost::phoenix::local_names;
      bp::function< phoenix::PHOENIX_MEMFUNC( dof ) > memfunc_dof;
      bp::function< phoenix::PHOENIX_MEMFUNC( dimensions ) > memfunc_dim;
      size_t result(0);
      boost::fusion::for_each
      ( 
        const_viewasref(*collapses_),
        bp::for_each
        ( 
          bpa::arg1,
          bp::lambda
          [
            let( bpl::_a = memfunc_dim( bpa::arg1 ) )
            [
              bp::if_( bpl::_a > bp::ref( dim ) )
                [ bp::ref(result) += memfunc_dof(bpa::arg1) ]
            ]
          ]
        ) 
      );
      return result;
    }
    INMANY( size_t ) :: dof() const 
    {
      namespace bp = boost::phoenix;
      namespace bpa = boost::phoenix::arg_names;
      bp::function< phoenix::PHOENIX_MEMFUNC( dof ) > memfunc;
      return boost::fusion::accumulate
             ( 
               const_viewasref(*collapses_), 0,
               bp::accumulate
               ( 
                 bpa::arg1, bpa::arg2,
                 bp::lambda[ memfunc(bpa::arg2) + bpa::arg1 ]
               ) 
             );
    }
    INMANY( void ) :: update( types::t_unsigned _d ) 
    {
      namespace bp = boost::phoenix;
      namespace bpa = boost::phoenix::arg_names;
      bp::function< phoenix::PHOENIX_MEMFUNC1( update ) > memfunc;
      boost::fusion::transform
      (
        viewasref(*collapses_),
        bp::for_each
        ( 
          bpa::arg1, 
          bp::lambda[ memfunc(bpa::arg1, bp::ref( _d ) ) ]
        ) 
      );
    }
    INMANY( void ) :: update_all() 
    {
      namespace bp = boost::phoenix;
      namespace bpa = boost::phoenix::arg_names;
      bp::function< phoenix::PHOENIX_MEMFUNC( update_all ) > memfunc;
      boost::fusion::transform
      (
        viewasref(*collapses_),
        bp::for_each
        ( 
          bpa::arg1, 
          bp::lambda[ memfunc(bpa::arg1) ]
        ) 
      );
    }
    INMANY( void ) :: reset() 
    {
      namespace bp = boost::phoenix;
      namespace bpa = boost::phoenix::arg_names;
      bp::function< phoenix::PHOENIX_MEMFUNC( reset ) > memfunc;
      boost::fusion::transform
      (
        viewasref(*collapses_),
        bp::for_each
        ( 
          bpa::arg1, 
          bp::lambda[ memfunc(bpa::arg1) ]
        ) 
      );
    }
    INMANY( void ) :: randomize( typename t_Vector :: value_type _howrandom ) 
    {
      namespace bp = boost::phoenix;
      namespace bpa = boost::phoenix::arg_names;
      bp::function< phoenix::PHOENIX_MEMFUNC1( randomize ) > memfunc;
      boost::fusion::transform
      (
        viewasref(*collapses_),
        bp::for_each
        ( 
          bpa::arg1, 
          bp::lambda[ memfunc(bpa::arg1, bp::val(_howrandom) ) ]
        ) 
      );
    }
    INMANY( size_t ) :: nbconfs() const 
    {
      namespace bp = boost::phoenix;
      namespace bpa = boost::phoenix::arg_names;
      bp::function< phoenix::PHOENIX_MEMFUNC( nbconfs ) > memfunc;
      return boost::fusion::accumulate
             ( 
               const_viewasref(*collapses_), 0,
               bp::accumulate
               ( 
                 bpa::arg1, bpa::arg2,
                 bp::lambda[ memfunc(bpa::arg2) + bpa::arg1 ]
               ) 
             );
    }
    INMANY( size_t ) :: dimensions() const 
    {
      namespace bp = boost::phoenix;
      namespace bpa = boost::phoenix::arg_names;
      namespace bpl = boost::phoenix::local_names;
      bp::function< phoenix::PHOENIX_MEMFUNC( dimensions ) > memfunc;
      size_t result(0);
      boost::fusion::for_each
      ( 
        const_viewasref(*collapses_),
        bp::for_each
        ( 
          bpa::arg1, 
          bp::lambda
          [ 
            bp::let( bpl::_a = memfunc( bpa::arg1 ) )
            [
              bp::if_( bpl::_a > bp::ref( result ) )
                [ bp::ref( result ) = bpl::_a ]
            ]
          ]
        ) 
      );
      return result;
    }

  template< class T_TRAITS >
  std::ostream& operator<<( std::ostream& _stream, const Many<T_TRAITS> &_many )
  {
    namespace bp = boost::phoenix;
    namespace bpa = boost::phoenix::arg_names;
    boost::fusion::for_each
    ( 
      _many.const_viewasref(*_many.collapses_),
      bp::for_each
      ( 
        bpa::arg1, 
        bp::lambda
        [
          bp::ref( _stream ) << bpa::arg1 << bp::val("\n")
        ]
      )
    );
    return _stream;
  }

  template< class T_STRUCTURES, class T_TRAITS >
   void init_many_collapses( const std::string &_desc, size_t _rank,
                             types::t_real _lambda, const T_STRUCTURES &_structures,
                             Many<T_TRAITS> &_many )
   {
     __DEBUGTRYBEGIN
     _many.mapping().init( _structures );
     boost::char_separator<char> sep(" ");
     typedef boost::tokenizer< boost::char_separator<char> > t_Tokenizer;
     t_Tokenizer notoken(_desc, sep);
     t_Tokenizer::const_iterator i_tok = notoken.begin();
     t_Tokenizer::const_iterator i_tok_end = notoken.end();
     const boost::regex rankrex("r=(\\d+)" );
     boost::match_results<std::string::const_iterator> what;
     size_t rank = _rank; 
     for(; i_tok != i_tok_end; ++i_tok )
     {
       CE::PosToConfs postoconfs( *Crystal::Structure::lattice );
       if( boost::regex_search( *i_tok, what, rankrex ) )
       {
         rank = boost::lexical_cast< size_t >( what.str(1) );
         continue;
       }
       __TRYCODE( postoconfs.create_positions( *i_tok );,
                  "Could not parse string " << _desc << "\n" )
       const size_t pos = _many.template addone<0>();
       _many.template collapses<0>( pos ).init( _structures, postoconfs );
       _many.template collapses<0>( pos ).regularization().lambda = _lambda;
       _many.template init<0>( rank, postoconfs.dof() );
     }
     __DEBUGTRYEND(,"Error while creating Many collapse/separables.\n" )
   }
  
  template< class T_MANY >
    void ManyState :: operator=( const T_MANY& _many )
    {
      coefficients_ = _many.coefficients();
      norms_.clear();
      boost::fusion::transform( *_many.collapses_, Save( norms_ ) );
    }
   template< class T_MANY > void ManyState :: reset( T_MANY& _many ) const
   {
     _many.coefficients() = coefficients_;
     boost::fusion::transform( *_many.collapses_, Reset( norms_ ) );
   }


} // end of CE namespace.

