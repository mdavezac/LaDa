#ifdef PYLADA_TYPEDEF
# error PYLADA_TYPEDEF already defined
#endif
#define PYLADA_TYPEDEF                                                          \
    (Pylada::math::rVector3d(*)( Pylada::math::rVector3d const&,                  \
                               Pylada::math::rMatrix3d const&,                  \
                               Pylada::math::rMatrix3d const& ))

#if PYLADA_CRYSTAL_MODULE != 1
  //! Refolds a periodic vector into the unit cell.
  PYLADA_INLINE math::rVector3d into_cell( math::rVector3d const &_vec, 
                                         math::rMatrix3d const &_cell, 
                                         math::rMatrix3d const &_inv)
    PYLADA_END(return (*PYLADA_TYPEDEF
                      api_capsule[PYLADA_SLOT(crystal)])(_vec, _cell, _inv);) 
#else
  api_capsule[PYLADA_SLOT(crystal)] = (void *)(PYLADA_TYPEDEF into_cell);
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(PYLADA_SLOT(crystal))
#include PYLADA_ASSIGN_SLOT(crystal)
  
#if PYLADA_CRYSTAL_MODULE != 1
  //! Refolds a periodic vector into the voronoi cell (eg first BZ or WignerSeitz).
  PYLADA_INLINE math::rVector3d into_voronoi( math::rVector3d const &_vec, 
                                            math::rMatrix3d const &_cell, 
                                            math::rMatrix3d const &_inv)
    PYLADA_END(return (*PYLADA_TYPEDEF
                      api_capsule[PYLADA_SLOT(crystal)])(_vec, _cell, _inv);) 
#else
  api_capsule[PYLADA_SLOT(crystal)] = (void *)(PYLADA_TYPEDEF into_voronoi);
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(PYLADA_SLOT(crystal))
#include PYLADA_ASSIGN_SLOT(crystal)

#if PYLADA_CRYSTAL_MODULE != 1
  //! \brief Refolds a periodic vector into a cell centered around zero (in
  //!        fractional coordinates).
  //! \details Since the vector is refolded in fractional coordinates, it may
  //!          or may not be the vector with smallest norm. Use math::rVector3d
  //!          into_voronoi() to get the equivalent vector with smallest norm.
  PYLADA_INLINE math::rVector3d zero_centered( math::rVector3d const &_vec, 
                                             math::rMatrix3d const &_cell, 
                                             math::rMatrix3d const &_inv)
    PYLADA_END(return (*PYLADA_TYPEDEF
                      api_capsule[PYLADA_SLOT(crystal)])(_vec, _cell, _inv);)
#else
  api_capsule[PYLADA_SLOT(crystal)] = (void *)(PYLADA_TYPEDEF zero_centered);
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(PYLADA_SLOT(crystal))
#include PYLADA_ASSIGN_SLOT(crystal)

#undef PYLADA_TYPEDEF

#if PYLADA_CRYSTAL_MODULE != 1
  //! Refolds a periodic vector into the unit cell.           
  inline math::rVector3d into_cell( math::rVector3d const &_vec, 
                                    math::rMatrix3d const &_cell )
    { return into_cell(_vec, _cell, _cell.inverse()); }
  //! \brief Refolds a periodic vector into the voronoi cell (eg first BZ or WignerSeitz).
  //! \details May fail if the matrix is a weird parameterization of the
  //!          lattice. It is best to use a grubber(...) cell . 
  inline math::rVector3d into_voronoi( math::rVector3d const &_vec, 
                                       math::rMatrix3d const &_cell )
    { return into_voronoi(_vec, _cell, _cell.inverse()); }

  //! \brief Refolds a periodic vector into a cell centered around zero (in
  //!        fractional coordinates).
  //! \details Since the vector is refolded in fractional coordinates, it may
  //!          or may not be the vector with smallest norm. Use math::rVector3d
  //!          into_voronoi() to get the equivalent vector with smallest norm.
  inline math::rVector3d zero_centered( math::rVector3d const &_vec, 
                                        math::rMatrix3d const &_cell )
    { return zero_centered(_vec, _cell, _cell.inverse()); }
#endif
