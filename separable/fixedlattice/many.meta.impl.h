//
//  Version: $Id$
//

//! \cond
namespace CE
{
// INCOLLAPSE( struct ) :: ApplyDof
// {
//   template< class T > size_t operator()( T &_t, size_t _r )
//   {
//     foreach( const typename T::value_type &_val, _t )
//       _r += _val.dof(); 
//     return _r;
//   }
// };
// INCOLLAPSE( struct ) :: ApplyNbconfs
// {
//   template< class T > size_t operator()( T &_t, size_t _r )
//   {
//     foreach( const typename T::value_type &_val, _t )
//       _r += _val.ApplyNbconfs(); 
//     return _r;
//   }
// };
  INCOLLAPSE( struct ) :: ApplyCurrentDof
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
  INCOLLAPSE2(template< class T_VECTOR > struct ) :: ApplyCreateAnB
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
      foreach( const typename T::value_type &_val, _t )
      {
        if( _val.dimensions() < self_.dim ) continue;
        const bblas::range colrange( rangestart,
                                     rangestart + _val.dof() );
    
        // create the X vector.
        bblas::vector_range< t_Vector > colX( X_, colrange );
        _val.create_X( i_, colX );
        rangestart += _v.dof();
      }
      rangestart = 0;
    }
  };
  INCOLLAPSE2(template< class T_MATRIX, class T_VECTOR > struct ) 
    :: ApplyRegularization
    {
      T_MATRIX &A_;
      T_VECTOR &b_;
      const Many& self_;
      size_t rangestart;
      ApplyRegularization   ( T_MATRIX &A_, T_VECTOR &_b, const Many& _self )
                          : A_(_A), b_(_b), self_(_self), rangestart(0) {}
      ApplyRegularization   ( const ApplyRegularization &_c )
                          : A_(_c.A_), b_(_c.b_), self_(_c.self_), rangestart(_c.rangestart) {}
      template< class T > void operator()( T &_t )
      {
        foreach( const typename T::value_type &_val, _t )
        {
          const bblas::range colrange( rangestart,
                                       rangestart + _val.dof() );
          bblas::matrix_range<T_MATRIX> sepmat( A_, seprange, seprange );
          bblas::vector_range<T_VECTOR> sepvec( b_, seprange );
          _val.regularization()( A_, b_, self_.dim ); 
      
          rangestart += _val.dof();
        }
      }
    };

  INCOLLAPSE( struct ) :: ApplyEvaluateOne
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

}
//! \endcond
