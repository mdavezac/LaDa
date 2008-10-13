//
//  Version: $Id$
//

namespace CE
{
# if defined( MANYHEAD ) || defined(INMANY) || defined(INMANY2)
#   error "Macros with same names."
# endif
# define MANYHEAD \
    ManyCollapses<T_TRAITS> 
#  define INMANY( var ) \
     template< class T_TRAITS >  var MANYHEAD
#  define INMANY2( code1, code2 ) \
     template< class T_TRAITS >  code1, code2 MANYHEAD

  INMANY2(template< class T_MATRIX, class T_VECTOR > void) 
    :: create_A_n_b( T_MATRIX &_A, T_VECTOR &_b )
    {
      namespace bblas = boost::numeric::ublas;
      __DEBUGTRYBEGIN
      // Loop over inequivalent configurations.
      const size_t cdof( current_dof() );
      t_Vector X( cdof );
      if( _A.size1() != cdof )  _A.resize( cdof, cdof );
      std::fill( _A.data().begin(), _A.data().end(), 
                 typename t_Matrix::value_type(0) );
      if( _b.size() != cdof )  _b.resize( cdof );
      std::fill( _b.data().begin(), _b.data().end(), 
                 typename t_Vector::value_type(0) );
      for( size_t i(0); i < mapping().size(); ++i )
      {
        // allows leave-one-out, or leave-many-out.
        if( mapping().do_skip(i) )  continue;
  
        std::fill( X.begin(), X.end(), typename t_Vector::value_type(0) );
        size_t range_start(0);
        foreach( t_FittingPair &pair, *fittingpairs_ )
        {
          t_Collapse &collapse( pair.first );
          if( collapse.dimensions() <= dim or collapse.dof() == 0 ) continue;
          const bblas::range colrange( rangestart,
                                       rangestart + collapse.dof() );
        
          // create the X vector.
          bblas::vector_range< t_Vector > colX( X, colrange );
          collapse.create_X( i, dim, colX );
          rangestart += collapse.dof();
        }
        
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
      create_A_n_b( _A, _b );

      // Apply regularization.
      size_t range_start(0);
      foreach( t_FittingPair &pair, *fittingpairs_ )
      {
        t_Collapse &collapse( pair.first );
        if( collapse.dimensions() <= dim ) continue;
        const bblas::range colrange( rangestart,
                                     rangestart + collapse.dof() );
        bblas::matrix_range<T_MATRIX> sepmat( _A, colrange, colrange );
        bblas::vector_range<T_VECTOR> sepvec( _b, colrange );
        collapse.regularization()( A_, _b, dim ); 
      
        rangestart += collapse.dof();
      }
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

  INMANY( typename MANYHEAD::t_Matrix::value_type ) :: evaluate( size_t _n) const 
  {
    types::t_real result(0);
    foreach( t_FittingPair &pair, *fittingpairs_ )
      result += pair.first.evaluate( _n ); 
    return result;
  }

  INMANY( size_t ) :: dof() const 
  {
    size_t result(0);
    foreach( t_FittingPair &pair, *fittingpairs_ )
      result += pair.first.dof();
    return result;
  }
  INMANY( size_t ) ::current_dof() const 
  {
    size_t result(0);
    foreach( t_FittingPair &pair, *fittingpairs_ )
      if( pair.first.dimensions() > dim ) 
        result += pair.first.dof();
    return result;
  }
  INMANY( void ) :: update( size_t _d ) const 
  {
    foreach( t_FittingPair &pair, *fittingpairs_ ) 
      if( pair.first.dimensions() > _d ) pair.update(_d);
  }
  INMANY( size_t ) :: nbconfs() const 
  {
    size_t result(0);
    foreach( t_FittingPair &pair, *fittingpairs_ )
      result += pair.first.nbconfs();
    return result;
  }
  INMANY( size_t ) :: dimensions() const 
  {
    size_t result(0);
    foreach( t_FittingPair &pair, *fittingpairs_ )
      result = std::max( pair.first.dimensions(), result );
    return result;
  }

  template< class T_TRAITS >
    std::ostream& operator<<( std::ostream& _stream,
                              const ManyCollapses<T_TRAITS> &_many )
    {
      foreach( t_FittingPair &pair, *fittingpairs_ )
        _stream << pair << "\n"; 
      return _stream;
    }

  INMANY( size_t ) :: addone()
  {
    const size_t N( fittingpairs_->size() );
    fittingpairs_->resize( N + 1 ); }
    return N;
  }
# ifdef MANYHEAD
#   undef MANYHEAD
# endif
# ifdef INMANY
#  undef INMANY
# endif
}
