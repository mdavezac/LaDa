//
//  Version: $Id$
//

//! \cond
#include <boost/concept_check.hpp>
#include <boost/concept/requires.hpp>
#include <boost/concept/usage.hpp>

namespace CE
{
  INMANY(template< class T_VECTOR > struct ) :: ApplyCreateAnB
  {
    typedef void result_type;
    T_VECTOR &X_;
    const size_t &i_;
    const size_t &dim;
    mutable size_t rangestart;
    ApplyCreateAnB   ( T_VECTOR &_X, const size_t _i, const size_t& _dim )
                   : X_(_X), i_(_i), dim(_dim), rangestart(0) {}
    ApplyCreateAnB   ( const ApplyCreateAnB & _c )
                   : X_(_c.X_), i_(_c.i_), dim(_c.dim),
                     rangestart(_c.rangestart) {}
    template< class T > result_type operator()( T &_t ) const
    {
      namespace bblas = boost::numeric::ublas;
      foreach( typename T::reference &_val, _t )
      {
        if( _val.dimensions() <= dim or _val.dof() == 0 ) continue;
        const bblas::range colrange( rangestart,
                                     rangestart + _val.dof() );
    
        // create the X vector.
        bblas::vector_range< t_Vector > colX( X_, colrange );
        _val.create_X( i_, dim, colX );
        rangestart += _val.dof();
      }
    }
  };
  INMANY2(template< class T_MATRIX, class T_VECTOR > struct ) 
    :: ApplyRegularization
    {
      typedef void result_type;
      T_MATRIX &A_;
      T_VECTOR &b_;
      const size_t& dim;
      mutable size_t rangestart;
      ApplyRegularization   ( T_MATRIX &_A, T_VECTOR &_b, const size_t& _dim )
                          : A_(_A), b_(_b), dim(_dim), rangestart(0) {}
      ApplyRegularization   ( const ApplyRegularization &_c )
                          : A_(_c.A_), b_(_c.b_), dim(_c.dim),
                            rangestart(_c.rangestart) {}
      template< class T > result_type operator()( T &_t ) const
      {
        namespace bblas = boost::numeric::ublas;
        foreach( typename T::const_reference &_val, _t )
        {
          if( _val.dimensions() <= dim ) continue;
          const bblas::range colrange( rangestart,
                                       rangestart + _val.dof() );
          bblas::matrix_range<T_MATRIX> sepmat( A_, colrange, colrange );
          bblas::vector_range<T_VECTOR> sepvec( b_, colrange );
          _val.regularization()( A_, b_, dim ); 
      
          rangestart += _val.dof();
        }
      }
    };

  INMANY( template< class T_COEFFICIENTS > struct )
    :: ApplyResize<const T_COEFFICIENTS>
  {
    typedef void result_type;
    mutable size_t rank_;
    T_COEFFICIENTS& coefficients_;
    ApplyResize   ( T_COEFFICIENTS & _coefficients ) 
                : coefficients_(_coefficients), rank_(0) {}
    ApplyResize   ( const ApplyResize & _c)
                : coefficients_(_c.coefficients_), rank_(_c.rank_) {}
    template< class T > result_type operator()( T &_t ) const
    {
      namespace bblas = boost::numeric::ublas;
      foreach( typename T::reference &_val, _t )
      {
        const size_t dof( _val.dof() );
        const bblas::range a( rank_, rank_ + dof );
        const bblas::range b( 0, _val.dimensions() );
        _val.coefficients_interface().set( coefficients_, a, b );
        rank_ += dof;
      }
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
    struct cadd_ref
    {
        template<typename Sig>
        struct result;
    
        template<typename U>
        struct result<cadd_ref(const U)>
        : boost::add_reference<const U>
        {};
    
        template <typename T>
        const T& operator()(const T& x) const
          { return x; }
    };
  }

  INMANY( template< class T > struct )
    ::  ViewAsRef : boost::fusion::transform_view 
                                   < T, details :: add_ref >
  {
    typedef boost::fusion::transform_view
            < T, details::add_ref > t_Base;

    explicit ViewAsRef( T& _v ) : t_Base( _v, details::add_ref() ) {}
  };
  INMANY( template< class T > struct )
    ::  cViewAsRef : boost::fusion::transform_view 
                                   < const T, details :: cadd_ref >
  {
    typedef boost::fusion::transform_view
            < const T, details::cadd_ref > t_Base;

    explicit cViewAsRef( const T& _v ) : t_Base( _v, details::cadd_ref() ) {}
  };

  struct ManyState :: Save
  {
    typedef void result_type;
    mutable ManyState :: t_Norms& norms;
    Save( ManyState::t_Norms& _norms ) : norms( _norms ) { norms.clear(); }
    Save( const Save &_c ) : norms( _c.norms ) {}
    template< class T >
      typename boost::enable_if
               <
                 ManyState::has_norms<typename T::const_reference>,
                 result_type 
               > :: type operator()( const T& _t ) const
        { 
          foreach( typename T::const_reference &_val, _t ) 
            norms.push_back( _val.norms );
        }
    template< class T >
      typename boost::disable_if
               <
                 ManyState::has_norms<typename T::const_reference>,
                 result_type 
               > :: type operator()( const T& _t ) const {}
  };

  struct ManyState :: Reset
  {
    typedef void result_type;
    mutable ManyState :: t_Norms :: const_iterator i_norm;
    ManyState :: t_Norms :: const_iterator i_norm_end;
    Reset   ( const ManyState::t_Norms& _norms ) 
          : i_norm( _norms.begin() ), i_norm_end( _norms.end() ) {}
    Reset( const Reset &_c ) : i_norm( _c.i_norm ), i_norm_end( i_norm_end ) {}


    template< class T >
      typename boost::enable_if
               <
                 ManyState::has_norms<typename T::reference>,
                 result_type 
               > :: type operator()( T& _t ) const
        { 
          foreach( typename T::reference &_val, _t ) 
          {
            __ASSERT( i_norm_end - i_norm <= 0, "Iterator out of range.\n" )
            _val.norms = *i_norm;
            ++i_norm;
          }
        }
    template< class T >
      typename boost::disable_if
               <
                 ManyState::has_norms<typename T::reference>,
                 result_type 
               > :: type operator()( T& _t ) const {}

  };

  template<class T > struct ManyState::has_norms
  {
    static const bool value = boost::remove_reference<T>::type::has_norms;
  };

}
//! \endcond
