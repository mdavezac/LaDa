#ifdef LADA_TYPEDEF
# error LADA_TYPEDEF already defined
#endif
#define LADA_TYPEDEF                                                          \
    (LaDa::math::rVector3d(*)( LaDa::math::rVector3d const&,                  \
                               LaDa::math::rMatrix3d const&,                  \
                               LaDa::math::rMatrix3d const& ))

#if LADA_CRYSTAL_MODULE != 1
  //! Refolds a periodic vector into the unit cell.
  LADA_INLINE math::rVector3d into_cell( math::rVector3d const &_vec, 
                                         math::rMatrix3d const &_cell, 
                                         math::rMatrix3d const &_inv)
    LADA_END( { return (*LADA_TYPEDEF
                         api_capsule[LADA_SLOT(crystal)])(_vec, _cell, _inv); } ) 
#else
  api_capsule[LADA_SLOT(crystal)] = (void *)(LADA_TYPEDEF into_cell);
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(crystal))
#include LADA_ASSIGN_SLOT(crystal)
  
#if LADA_CRYSTAL_MODULE != 1
  //! Refolds a periodic vector into the voronoi cell (eg first BZ or WignerSeitz).
  LADA_INLINE math::rVector3d into_voronoi( math::rVector3d const &_vec, 
                                            math::rMatrix3d const &_cell, 
                                            math::rMatrix3d const &_inv)
    LADA_END({ return (*LADA_TYPEDEF
                        api_capsule[LADA_SLOT(crystal)])(_vec, _cell, _inv); }) 
#else
  api_capsule[LADA_SLOT(crystal)] = (void *)(LADA_TYPEDEF into_voronoi);
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(crystal))
#include LADA_ASSIGN_SLOT(crystal)

#if LADA_CRYSTAL_MODULE != 1
  //! \brief Refolds a periodic vector into a cell centered around zero (in
  //!        fractional coordinates).
  //! \details Since the vector is refolded in fractional coordinates, it may
  //!          or may not be the vector with smallest norm. Use math::rVector3d
  //!          into_voronoi() to get the equivalent vector with smallest norm.
  LADA_INLINE math::rVector3d zero_centered( math::rVector3d const &_vec, 
                                             math::rMatrix3d const &_cell, 
                                             math::rMatrix3d const &_inv)
    LADA_END({ return (*LADA_TYPEDEF
                        api_capsule[LADA_SLOT(crystal)])(_vec, _cell, _inv); })
#else
  api_capsule[LADA_SLOT(crystal)] = (void *)(LADA_TYPEDEF zero_centered);
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(crystal))
#include LADA_ASSIGN_SLOT(crystal)

#undef LADA_TYPEDEF

#if LADA_CRYSTAL_MODULE != 1
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
