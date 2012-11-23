#if LADA_CRYSTAL_MODULE != 1

namespace LaDa
{
  namespace crystal 
  {
    namespace 
    {
#endif

#if LADA_CRYSTAL_MODULE != 1
  //! \brief Returns true if two structures are equivalent. 
  //! \details Two structures are equivalent in a crystallographic sense,
  //!          e.g. without reference to cartesian coordinates or possible
  //!          motif rotations which leave the lattice itself invariant. A
  //!          supercell is *not* equivalent to its lattice, unless it is a
  //!          trivial supercell.
  //! \param[in] _a: The first structure.
  //! \param[in] _b: The second structure.
  //! \param[in] scale: whether to take the scale into account. Defaults to true.
  //! \param[in] cartesian: whether to take into account differences in
  //!            cartesian coordinates. Defaults to true. If False, then
  //!            comparison is according to mathematical definition of a
  //!            lattice. If True, comparison is according to
  //!            crystallographic comparison.
  //! \param[in] tolerance: Tolerance when comparing distances. Defaults to
  //!            types::t_real. It is in the same units as the structures scales, if
  //!            that is taken into account, otherwise, it is in the same
  //!            units as _a.scale.
  bool equivalent( Structure const &_a, Structure const &_b,
                   bool with_scale=true, bool with_cartesian=true,
                   types::t_real _tol = types::tolerance )
    LADA_END({ return (*(bool(*)(Structure const&, Structure const &, bool, bool, types::t_real))
                       api_capsule[LADA_SLOT(crystal)])(_a, _b, with_scale, with_cartesian, _tol); })
#else
  api_capsule[LADA_SLOT(crystal)] = (void *)equivalent;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(crystal))
#include LADA_ASSIGN_SLOT(crystal)

#if LADA_CRYSTAL_MODULE != 1
    } 
  } // namespace crystal
} // namespace LaDa
#endif
