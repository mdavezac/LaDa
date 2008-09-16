//
//  Version: $Id$
//

//! \cond
#include <boost/lambda/lambda.hpp>

namespace CE
{
  INMANY( struct ) :: ApplyDimensions
  {
    typedef size_t result_type;
    template< class T > result_type operator()( const T &_t, size_t _r )
    {
      foreach( const typename boost::remove_pointer<typename T::value_type> :: type &_val, _t )
        _r = std::max( _r, _val.dimensions() );
      return _r;
    }
  };
  INMANY( struct ) :: ApplyCurrentDof
  {
    typedef size_t result_type;
    const size_t &dim_;
    ApplyCurrentDof( const size_t &_dim ) : dim_(_dim) {}
    ApplyCurrentDof( const ApplyCurrentDof &_c ) : dim_(_c.dim_) {}
    template< class T > result_type operator()( T &_t, size_t _r )
    {
      foreach( const typename T::value_type &_val, _t )
        if( _val.dimensions() < dim_ ) _r += _val.dof(); 
      return _r;
    }
  };
  INMANY(template< class T_VECTOR > struct ) :: ApplyCreateAnB
  {
    typedef void result_type;
    T_VECTOR &X_;
    const size_t &i_;
    const Many& self_;
    size_t rangestart;
    ApplyCreateAnB   ( T_VECTOR &_X, const size_t _i, const Many& _self )
                   : X_(_X), i_(_i), self_(_self), rangestart(0) {}
    ApplyCreateAnB   ( const ApplyCreateAnB & _c )
                   : X_(_c.X_), i_(_c.i_), self_(_c.self_), rangestart(_c.rangestart) {}
    template< class T > result_type operator()( T &_t )
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
      typedef void result_type;
      T_MATRIX &A_;
      T_VECTOR &b_;
      const Many& self_;
      size_t rangestart;
      ApplyRegularization   ( T_MATRIX &_A, T_VECTOR &_b, const Many& _self )
                          : A_(_A), b_(_b), self_(_self), rangestart(0) {}
      ApplyRegularization   ( const ApplyRegularization &_c )
                          : A_(_c.A_), b_(_c.b_), self_(_c.self_),
                            rangestart(_c.rangestart) {}
      template< class T > result_type operator()( T &_t ) 
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
    typedef typename t_Matrix::value_type result_type;
    size_t &n_;
    ApplyEvaluateOne( size_t &_n ) : n_(_n) {}
    ApplyEvaluateOne( const ApplyEvaluateOne &_c ) : n_(_c.n_) {}
    template< class T > result_type operator()( T &_t, result_type _r )
    {
      foreach( const typename T::value_type &_val, _t )
        _r += _val.evaluate(n_);
      return _r;
    }
  };

  INMANY( template< class T_COEFFICIENTS > struct )
    :: ApplyResize<const T_COEFFICIENTS>
  {
    typedef void result_type;
    size_t rank_;
    const T_COEFFICIENTS& coefficients_;
    ApplyResize   ( const T_COEFFICIENTS & _coefficients ) 
                : coefficients_(_coefficients), rank_(0) {}
    ApplyResize   ( const ApplyResize & _c)
                : coefficients_(_c.coefficients_), rank_(_c.rank_) {}
    template< class T > result_type operator()( T &_t ) const
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
    typedef T_STREAM& result_type;
    mutable T_STREAM &stream_;
    PrintToStream( T_STREAM &_stream ) : stream_(_stream) {}
    PrintToStream( PrintToStream &_c ) : stream_(_c.stream_) {}
    template< class T > result_type operator()( const T &_t ) const
    {
      namespace bl = boost::lambda;
      foreach( const typename T::value_type &_val, _t )
        stream_  << _val << "\n";
      return stream_;
    }
  };

  namespace details
  {
    struct add_ref
    {
        template<typename Sig>
        struct result;
    
        template<typename U>
        struct result<add_ref(U)>
        : boost::add_reference<U>
        {};
    
        template <typename T>
        T& operator()(T& x) const
          { return x; }
    };
  }

  INMANY( template< class T > struct ) ::  ViewAsRef : boost::fusion::transform_view
              < 
                T, details :: add_ref
              >
  {
    typedef boost::fusion::transform_view
            < T, details::add_ref > t_Base;

    explicit ViewAsRef( T& _v ) : t_Base( _v, details::add_ref() ) {}
  };

  struct ManyState :: Save
  {
    typedef void result_type;
    ManyState :: t_Norms& norms;
    Save( ManyState::t_Norms& _norms ) : norms( _norms ) { norms.clear(); }
    Save( const Save &_c ) : norms( _c.norms ) {}
    template< class T >
      result_type operator()( const T& _t )
        { foreach( const typename T::value_type &_val, _t ) call( _t ); }
    template< class T >
      typename boost::enable_if< ManyState::has_norms<T>, result_type > :: type
        call( const T& _t ) { norms.push_back( _t.norms ); }
    template< class T >
      typename boost::disable_if< ManyState::has_norms<T>, result_type > :: type
        call( const T& _t ) {}
  };

  struct ManyState :: Reset
  {
    typedef void result_type;
    ManyState :: t_Norms :: const_iterator i_norm;
    ManyState :: t_Norms :: const_iterator i_norm_end;
    Reset   ( const ManyState::t_Norms& _norms ) 
          : i_norm( _norms.begin() ), i_norm_end( _norms.end() ) {}
    Reset( const Reset &_c ) : i_norm( _c.i_norm ), i_norm_end( i_norm_end ) {}
    template< class T >
      result_type operator()( T& _t )
          { foreach( typename T::value_type &_val, _t ) call( _t ); }
    template< class T >
      typename boost::enable_if< ManyState::has_norms<T>, result_type > :: type
        call( T& _t )
        {
          __ASSERT( i_norm_end - i_norm > 0, "Iterator out of range.\n" )
          _t.norms = *i_norm;
          ++i_norm;
        }
    template< class T >
      typename  boost::disable_if< ManyState::has_norms<T>, result_type > :: type
        call( T& _t ) {}
  };

  template<class T > class ManyState::has_norms
  {
    template<class TT> struct CheckForNorms
    {
      BOOST_CONCEPT_USAGE(CheckForNorms)
      {
        t.norms.clear();
        t.norms.begin();
        t.norms.end();
      }
      TT t;
    };
    typedef char One;
    struct Two { One two[2]; };
    template<class U> One (&test( U*, CheckForNorms<U>* a=0 ));
    Two (&test(...));
    enum { value = 1 == sizeof( has_norms::test( (T*)0 ) ) };
  };


}
//! \endcond
