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

#include "prepare.h"

namespace CE
{
# if defined( COLHEAD ) || defined(INCOLLAPSE) || defined(INCOLLAPSE2)
#   error "Macros with same names."
# endif
# define COLHEAD \
    Many<T_TRAITS> 
#  define INCOLLAPSE( var ) \
     template< class T_TRAITS >  var COLHEAD
#  define INCOLLAPSE2( code1, code2 ) \
     template< class T_TRAITS >  code1, code2 COLHEAD



    INCOLLAPSE2(template< class T_MATRIX, class T_VECTOR > void) 
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
        const t_Collapse *first_col = collapses_->front();
        for( size_t i(0); i < mapping().size(); ++i )
        {
          // allows leave-one-out, or leave-many-out.
          if( mapping().do_skip(i) )  continue;

          bf::for_each( *collapses_, ApplyCreateAnB<t_Vector>( X, i, *this ) );
          
          _A += mapping().weight(i) * bblas::outer_prod( X, X ); 
          _b += mapping().weight(i) * mapping().target(i) * X;
        }
        __DEBUGTRYEND(, "Error in Many::create_A_n_b()\n" )
      }

    INCOLLAPSE2( template< class T_MATRIX, class T_VECTOR > void )
      :: operator()( T_MATRIX &_A, T_VECTOR &_b, types::t_unsigned _dim )
      {
        __DEBUGTRYBEGIN
        namespace bblas = boost::numeric::ublas;
        dim = _dim;
        foreach( t_Collapse &collapse, *collapses_ ) collapse.self()->dim = _dim;
        create_A_n_b( _A, _b );

        ApplyRegularization< T_MATRIX, T_VECTOR > applyreg( _A, _b, *this );
        boost::fusion::for_each( *collapses_, applyreg )
        __DEBUGTRYEND(, "Error in Many::operator()()\n" )
      }
    
    INCOLLAPSE( opt::ErrorTuple ) :: evaluate() const 
    {
      __DEBUGTRYBEGIN
      opt::ErrorTuple error;
      const t_Collapse &first_collapse = collapses_->front();
      for(size_t n(0); n < mapping().size(); ++n )
        if( not mapping().do_skip(n) )
          error += opt::ErrorTuple( mapping().target( n ) - evaluate(n),
                                    mapping().weight( n ) );
      return error;
      __DEBUGTRYEND(, "Error in Many::evaluate()\n" )
    }

    INCOLLAPSE2( template< class T_COLLAPSE, class T_SEPARABLES > size_t )
      ::  add_as_is()
      {
        namespace bm = boost::mpl;
        namespace bg = boost::fusion;
        typedef make_ptrlist::apply<T_COLLAPSE> :: type t_ColPtrList;
        typedef make_ptrlist::apply<T_SEPARABLES> :: type t_SepPtrList;
        bf::filter_view< t_ListOfSeparables,
                         bm::is_same<_, t_ColPtrList > > colview( *collapse_ );
        bf::filter_view< t_ListOfSeparables,
                         bm::is_same<_, t_SepPtrList > > sepview( *separables_ );
        bf::at<0>( colview ).push_back( new T_COLLAPSE );
        bf::at<0>( sepview ).push_back( new T_SEPARABLES );
        bf::at<0>( colview ).back().init( &bf::at<0>( sepview ).back() );
        return bf::at<0>( sepview ).size() - 1;
      }
    INCOLLAPSE2( template< class T_COLLAPSE, class T_SEPARABLES > size_t )
      ::  wrap_n_add()
      {
        typedef typename SeparablesWithMatrixRange<t_OrigSeparables>
                                                  :: other t_newSeparables;
        typedef typename CollapseWithNewSeparables<t_Separables> 
                                                  :: other t_newCollapse;
        return add_as_is< t_newSeparables, t_newCollapse >();
      }
    INCOLLAPSE( void ) :: init( size_t _index, size_t _rank, size_t _dimensions )
    {
      namespace bm = boost::numeric::ublas;
      size_t rank( _rank );
      size_t dimensions( _dimensions );
      foreach( const t_Collapse &collapse, *collapses_ )
      {
        rank = std::max( rank, collapse.self().dof() );
        dimensions = std::max( dimensions, collapse.self().dimensions() );
      }
      coefficients_.resize( rank, dimensions );
      rank = 0; dimensions = 0;
      foreach( const t_Collapse &collapse, *collapses_ )
      {
        const thisdof( collapse.self().dof() );
        const thisdim( collapse.self().dimensions() );
        const bm a( rank, rank + thisdof );
        const bm b( dimensions, dimensions + thisdim );
        collapse.self().coefficients_interface().set( coefficients, a, b );
        rank += thisdof;
        dimensions += thisdim;
      }
    }

#  undef COLHEAD
#  undef INCOLLAPSE
#  undef INCOLLAPSE2

  template< class T_TRAITS >
  std::ostream& operator<<( std::ostream& _stream, const Many<T_TRAITS> &_many )
  {
    for( size_t i(0); i < _many.size(); ++i )
      _stream << _many.separables(i) << "\n";
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
       const size_t index = _many.addone<T_COLLAPSE, T_SEPARABLES>();
       _many.collapse( index )->init( _structures, postoconfs );
       _many.collapse( index )->regularization().lambda = _lambda;
       _many.init( index, rank, postoconfs.dof() );
     }
     __DEBUGTRYEND(,"Error while creating Many collapse/separables.\n" )
   }
} // end of CE namespace.

