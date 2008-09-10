//
//  Version: $Id$
//
#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/matrix.hpp>
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/numeric/ublas/operation.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

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

    INCOLLAPSE( struct ) :: make_ptrvector
    {
      typedef 

    };
    INCOLLAPSE( size_t ) :: dimensions() const 
    {
      size_t result(0);
      foreach( const t_Collapse *collapse, *collapses_ )
        result += collapse->self_()->self()->dimensions();
      return result;
    } 
    INCOLLAPSE( size_t ) :: dof() const 
    {
      size_t result(0);
      foreach( const t_Collapse &collapse, *collapses_ )
        result += collapse.self()->dof();
      return result;
    }
    INCOLLAPSE( size_t ) :: nbconfs() const 
    {
      size_t result(0);
      foreach( const t_Collapse &collapse, *collapses_ )
        result += collapse.self()->nbconfs();
      return result;
    } 
    INCOLLAPSE( size_t ) :: current_dof() const 
    {
      size_t result(0);
      foreach( const t_Collapse &collapse, *collapses_ )
        if( collapse.self()->dimensions() < dim ) result += collapse.self()->dof();
    }
    INCOLLAPSE( void ) :: reset() 
    {
      foreach( const t_Collapse &collapse, *collapses_ )
        collapse.self()->reset(); 
    }

    INCOLLAPSE2(template< class T_MATRIX, class T_VECTOR > void) 
      :: create_A_n_b( T_MATRIX &_A, T_VECTOR &_b )
      {
        __DEBUGTRYBEGIN
        namespace bblas = boost::numeric::ublas;
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

          size_t rangestart(0);
          foreach( t_Collapse &_collapse, *collapses_ )
          {
            if( first_col->self()->dimensions() < dim ) continue;
            const bblas::range colrange( rangestart,
                                         rangestart + first_col->self()->dof() );
      
            // create the X vector.
            bblas::vector_range< t_Vector > colX( X, colrange );
            std::fill( colX.begin(), colX.end(), typename t_Vector::value_type(0) );
            collapse.create_X( i, colX );
            rangestart += first_col->self()->dof();
          }
          
      
          _A += first_col->self()->weight(i) * bblas::outer_prod( X, X ); 
          _b += first_col->self()->weight(i) * mapping().target(i) * X;
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

        size_t rangestart(0);
        foreach( t_Collapse &_collapse, *collapses_ )
        {
          const bblas::range seprange( rangestart,
                                       rangestart + collapse.self()->dof() );
          bblas::matrix_range<T_MATRIX> sepmat( _A, seprange, seprange );
          bblas::vector_range<T_VECTOR> sepvec( _b, seprange );
          collapse.self()->regularization()( _A, _b, _dim); 
          rangestart += collapse.self()->dof();
        }
        __DEBUGTRYEND(, "Error in Many::operator()()\n" )
      }
    
    INCOLLAPSE( void ) :: randomize( typename t_Vector :: value_type _howrandom )
    {
      foreach( t_Collapse &collapse, *collapses_ )
        collapse.self()->randomize( _howrandom );
    }


    INCOLLAPSE( void ) :: update_all()
    {
      __DEBUGTRYBEGIN
      foreach( t_Collapse &collapse, *collapses_ ) collapse.self()->update_all();
      __DEBUGTRYEND(, "Error in Many::update_all()\n" )
    }
    INCOLLAPSE( void ) :: update( types::t_unsigned _d )
    {
      __DEBUGTRYBEGIN
      foreach( t_Collapse &collapse, *collapses_ ) collapse.self()->update(_d);
      __DEBUGTRYEND(, "Error in Many::update()\n" )
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

    INCOLLAPSE( typename COLHEAD::t_Matrix::value_type ) :: evaluate(size_t _n) const 
    {
      __DEBUGTRYBEGIN
      typename t_Matrix :: value_type intermed(0);
      foreach( t_Collapse &collapse, *collapses_ )
        intermed += collapse.self()->evaluate(_n);

      return intermed;
      __DEBUGTRYEND(, "Error in Many::evaluate( size_t )\n" )
    }

    INCOLLAPSE2( template< class T_COLLAPSE, class T_SEPARABLES > size_t )
      ::  add_as_is()
      {
        separables_->push_back( new opt::Indirection< T_SEPARABLES > );
        collapses_->push_back( new opt::Indirection< T_COLLAPSE > );
        collapses_->back().self().init( separables_->back().self() );
        return separables_->size() - 1;
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

