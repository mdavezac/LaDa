//
//  Version: $Id$
//

//! \cond

namespace CE
{
  INMANY( struct ) :: ApplyDimensions
  {
    template< class T > size_t operator()( const T &_t, size_t _r )
    {
      foreach( const typename T::value_type &_val, _t )
        _r = std::max( _r, _val.dimensions() );
      return _r;
    }
  };
  INMANY( struct ) :: ApplyCurrentDof
  {
    size_t &dim_;
    ApplyCurrentDof( size_t &_dim ) : dim_(_dim) {}
    ApplyCurrentDof( const ApplyCurrentDof &_c ) : dim_(_c.dim_) {}
    template< class T > size_t operator()( T &_t, size_t _r )
    {
      foreach( const typename T::value_type &_val, _t )
        if( _val.dimensions() < dim_ ) _r += _val.dof(); 
      return _r;
    }
  };
  INMANY(template< class T_VECTOR > struct ) :: ApplyCreateAnB
  {
    T_VECTOR &X_;
    const size_t &i_;
    const Many& self_;
    size_t rangestart;
    ApplyCreateAnB   ( T_VECTOR &_X, const size_t _i, const Many& _self )
                   : X_(_X), i_(_i), self_(_self), rangestart(0) {}
    ApplyCreateAnB   ( const ApplyCreateAnB & _c )
                   : X_(_c.X_), i_(_c.i_), self_(_c.self), rangestart(_c.rangestart) {}
    template< class T > void operator()( T &_t )
    {
      namespace bblas = boost::numeric::ublas;
      foreach( const typename T::value_type &_val, _t )
      {
        if( _val.dimensions() < self_.dim ) continue;
        const bblas::range colrange( rangestart,
                                     rangestart + _val.dof() );
    
        // create the X vector.
        bblas::vector_range< t_Vector > colX( X_, colrange );
        _val.create_X( i_, colX );
        rangestart += _val.dof();
      }
      rangestart = 0;
    }
  };
  INMANY2(template< class T_MATRIX, class T_VECTOR > struct ) 
    :: ApplyRegularization
    {
      T_MATRIX &A_;
      T_VECTOR &b_;
      const Many& self_;
      size_t rangestart;
      ApplyRegularization   ( T_MATRIX &_A, T_VECTOR &_b, const Many& _self )
                          : A_(_A), b_(_b), self_(_self), rangestart(0) {}
      ApplyRegularization   ( const ApplyRegularization &_c )
                          : A_(_c.A_), b_(_c.b_), self_(_c.self_), rangestart(_c.rangestart) {}
      template< class T > void operator()( T &_t )
      {
        namespace bblas = boost::numeric::ublas;
        foreach( const typename T::value_type &_val, _t )
        {
          const bblas::range colrange( rangestart,
                                       rangestart + _val.dof() );
          bblas::matrix_range<T_MATRIX> sepmat( A_, colrange, colrange );
          bblas::vector_range<T_VECTOR> sepvec( b_, colrange );
          _val.regularization()( A_, b_, self_.dim ); 
      
          rangestart += _val.dof();
        }
      }
    };

  INMANY( struct ) :: ApplyEvaluateOne
  {
    size_t &n_;
    ApplyEvaluateOne( size_t &_n ) : n_(_n) {}
    ApplyEvaluateOne( const ApplyEvaluateOne &_c ) : n_(_c.n_) {}
    template< class T > size_t operator()( T &_t, typename t_Matrix::value_type _r )
    {
      foreach( const typename T::value_type &_val, _t )
        _r += _val.evaluate(n_);
      return _r;
    }
  };

  INMANY( template< class T_COEFFICIENTS > struct ) :: ApplyResize
  {
    size_t rank_;
    T_COEFFICIENTS& coefficients_;
    ApplyResize( T_COEFFICIENTS & _coefficients ) : coefficients_(_coefficients), rank_(0) {}
    ApplyResize( const ApplyResize & _c) : coefficients_(_c/coefficients_), rank_(_c.rank_) {}
    template< class T > size_t operator()( T &_t, typename t_Matrix::value_type _r )
    {
      namespace bblas = boost::numeric::ublas;
      foreach( typename T::value_type &_val, _t )
      {
        const size_t dof( _val.dof() );
        const bblas::range a( rank_, rank_ + dof );
        const bblas::range b( 0, _val.dimensions() );
        _val.coefficients_interface().set( coefficients_, a, b );
        rank_ += dof;
      }
    }
  };

  INMANY( template< class T_STREAM > struct ) :: PrintToStream
  {
    T_STREAM &stream_;
    PrintToStream( T_STREAM &_stream ) : stream_(_stream) {}
    PrintToStream( PrintToStream &_c ) : stream_(_c.stream_) {}
    template< class T > void operator()( const T &_t )
    {
      foreach( typename T::value_type &_val, _t )
        stream_ << _val << "\n";
      return stream_;
    }
  };
}
//! \endcond
