#if PYLADA_MATH_MODULE != 1

  //! Exponentiation for integer powers.
  template<class T>
    inline T ipow(T base, size_t exponent) 
    {
      T result(1);
      for(; exponent; --exponent, result *= base);
      return result;
    }

# ifdef PYLADA_MACRO
#   error PYLADA_MACRO already defined.
# endif
# define PYLADA_MACRO(name, default_) \
    template<class T_DERIVED>                                                   \
      inline bool name( Eigen::DenseBase<T_DERIVED> const & _in,                \
                        typename Eigen::DenseBase<T_DERIVED>::RealScalar const &_tol ) \
      {                                                                         \
        for(int i(0); i < _in.rows(); ++i)                                   \
          for(int j(0); j < _in.cols(); ++j)                                 \
            if(not name(_in.coeff(i,j), _tol)) return false;                    \
        return true;                                                            \
      }                                                                         \
    template<class T_DERIVED>                                                   \
      inline bool name(Eigen::DenseBase<T_DERIVED> const & _in)                 \
      {                                                                         \
        for(int i(0); i < _in.rows(); ++i)                                   \
          for(int j(0); j < _in.cols(); ++j)                                 \
            if(not name(_in.coeff(i,j))) return false;                          \
        return true;                                                            \
      }
  PYLADA_MACRO(is_integer, types::t_real _tol = types::tolerance)
  PYLADA_MACRO(is_null,    types::t_real _tol = types::tolerance)
# undef PYLADA_MACRO
  //! True if an eigen array or matrix is the identity.
  template<class T_DERIVED>
    inline bool is_identity( Eigen::DenseBase<T_DERIVED> const & _in,
                             typename Eigen::DenseBase<T_DERIVED>::RealScalar const &_tol )
    {
      if(_in.cols() == 1)
      {
        for(int j(0); j < _in.cols(); ++j)
          if(not is_identity(_in.coeff(j,0), _tol)) return false;
      }
      else 
      {
        for(int i(0); i < _in.rows(); ++i)
          for(int j(0); j < _in.cols(); ++j)
            if(i==j) { if( not is_identity(_in.coeff(i,j), _tol)) return false; }
            else if(not is_null(_in.coeff(i,j), _tol)) return false;
      }
      return true;
    }
  //! True if an eigen array or matrix is the identity.
  template<class T_DERIVED>
    inline bool is_identity(Eigen::DenseBase<T_DERIVED> const & _in)
    {
      if(_in.cols() == 1)
      {
        for(int j(0); j < _in.rows(); ++j)
          if(not is_identity(_in.coeff(j,0))) return false;
      }
      else 
      {
        for(int i(0); i < _in.rows(); ++i)
          for(int j(0); j < _in.cols(); ++j)
            if(i==j) { if( not is_identity(_in.coeff(i,j))) return false; }
            else if(not is_null(_in.coeff(i,j))) return false;
      }
      return true;
    }
  //! True if two eigen arrays or matrices are equal, according to math::eq(). 
  template<class T_DERIVED0, class T_DERIVED1> 
    inline bool eq( Eigen::DenseBase<T_DERIVED0> const & _a,
                    Eigen::DenseBase<T_DERIVED1> const & _b,
                    typename Eigen::DenseBase<T_DERIVED0>::RealScalar const &_tol )
    {
      typedef typename Eigen::DenseBase<T_DERIVED0>::RealScalar t_real0;
      typedef typename Eigen::DenseBase<T_DERIVED1>::RealScalar t_real1;
      //BOOST_STATIC_ASSERT((boost::is_same<t_real0, t_real1>::value));
#     ifdef PYLADA_DEBUG
        if(_a.rows() != _b.rows() or _a.cols() != _b.cols())
          BOOST_THROW_EXCEPTION(error::array_of_different_sizes());
#     endif
      for(int i(0); i < _a.rows(); ++i)
        for(int j(0); j < _a.cols(); ++j)
          if(not eq(_a.coeff(i,j), _b.coeff(i,j), _tol)) return false;
      return true;
    }
  //! True if two eigen arrays or matrices are equal, according to math::eq(). 
  template<class T_DERIVED0, class T_DERIVED1>
    inline bool eq( Eigen::DenseBase<T_DERIVED0> const & _a,
                    Eigen::DenseBase<T_DERIVED1> const & _b )
    {
      typedef typename Eigen::DenseBase<T_DERIVED0>::RealScalar t_real0;
      typedef typename Eigen::DenseBase<T_DERIVED1>::RealScalar t_real1;
      //BOOST_STATIC_ASSERT((boost::is_same<t_real0, t_real1>::value));
#     ifdef PYLADA_DEBUG
        if(_a.rows() != _b.rows() or _a.cols() != _b.cols())
          BOOST_THROW_EXCEPTION(error::array_of_different_sizes());
#     endif
      for(int i(0); i < _a.rows(); ++i)
        for(int j(0); j < _a.cols(); ++j)
          if(not eq(_a.coeff(i,j), _b.coeff(i,j))) return false;
      return true;
    }
  //! True if two symmetry operations are equal.
  inline bool eq( Affine3d const & _a, Affine3d const & _b )
    { return eq(_a.linear(), _b.linear()) and eq(_a.translation(), _b.translation()); }
  //! True if two symmetry operations are equal.
  inline bool eq( Affine3d const & _a, Affine3d const & _b, Affine3d::Scalar const &_tol )
    { return eq(_a.linear(), _b.linear(), _tol) and eq(_a.translation(), _b.translation(), _tol); }
  //! True if two eigen arrays or matrices are not equal, according to math::eq(). 
  template<class T_DERIVED0, class T_DERIVED1> 
    inline bool neq( Eigen::DenseBase<T_DERIVED0> const & _a,
                     Eigen::DenseBase<T_DERIVED1> const & _b,
                     typename Eigen::DenseBase<T_DERIVED0>::RealScalar const &_tol )
      { return not eq(_a, _b, _tol); }
  //! True if two eigen arrays or matrices are not equal, according to math::eq(). 
  template<class T_DERIVED0, class T_DERIVED1> 
    inline bool neq( Eigen::DenseBase<T_DERIVED0> const & _a,
                     Eigen::DenseBase<T_DERIVED1> const & _b )
      { return not eq(_a, _b); }
  //! True if two symmetry operations are not equal.
  inline bool neq( Affine3d const & _a, Affine3d const & _b ) { return not eq(_a, _b); }
  //! True if two symmetry operations are not equal.
  inline bool neq( Affine3d const & _a, Affine3d const & _b, Affine3d::Scalar const &_tol )
      { return not eq(_a, _b, _tol); }

  //! True if two vectors are periodic images with respect to a cell.
  inline bool are_periodic_images( math::rVector3d const &_a,
                                   math::rVector3d const &_b, 
                                   math::rMatrix3d const &_inv_cell )
   { return math::is_integer(_inv_cell * (_a - _b)); }
  //! True if two vectors are periodic images with respect to a cell.
  inline bool are_periodic_images( math::rVector3d const &_a,
                                   math::rVector3d const &_b, 
                                   math::rMatrix3d const &_inv_cell,
                                   types::t_real _tol )
   { return math::is_integer(_inv_cell * (_a - _b), _tol); }


  //! Defines a floor function.
  inline rMatrix3d floor(rMatrix3d const &_matrix)
  {
    rMatrix3d result;
    for(int i(0); i < 3; ++i)
      for(int j(0); j < 3; ++j)
        result(i,j) = std::floor(_matrix(i,j));
    return result;
  }
  //! Defines a floor function.
  inline math::rVector3d floor(math::rVector3d const &_v) 
  {
    return math::rVector3d( std::floor(_v(0)), 
                            std::floor(_v(1)), 
                            std::floor(_v(2)) );
  }

  //! Rounding off.
  inline types::t_real round(types::t_real r) 
    { return (r > 0.0) ? std::floor(r + 0.5) : std::ceil(r - 0.5); }

  //! Rounding off.
  inline math::rVector3d round(math::rVector3d const &_v)
    { return math::rVector3d( round(_v(0)), round(_v(1)), round(_v(2)) ); }

  //! Maximum component.
  template<class T> inline T max( Eigen::Matrix<T, 3, 1>const &_v)
    { return std::max(std::max(_v(0), _v(1)), _v(2)); }

  //! Casts to lower integer accounting for numerical noise.
  inline iVector3d floor_int(rVector3d const &_t )
    { return iVector3d(floor_int(_t(0)), floor_int(_t(1)), floor_int(_t(2))); }
  
  //! Absolute norm.
  inline rVector3d::Scalar absnorm(rVector3d const &_a)
    { return std::max(std::max(std::abs(_a[0]), std::abs(_a[1])), std::abs(_a[2])); }

#endif
