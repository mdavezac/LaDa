#if LADA_CRYSTAL_MODULE != 1
namespace LaDa
{
  namespace crystal
  {
    namespace
    {
#endif

#if LADA_CRYSTAL_MODULE != 1
  //! \brief Creates list of coordination shells up to given order.
  //! \returns A list of lists of tuples. The outer list is over coordination shells.
  //!          The inner list references the atoms in a shell.
  //!          Each innermost tuple contains a reference to the atom in question,
  //!          a translation vector to its periodic image inside the relevant shell, 
  //!          and the distance from the center to the relevant periodic image.
  //! \param[in] _structure : Structure for which to determine coordination shells.
  //! \param[in] _nshells : Number of shells to compute.
  //! \param[in] _center : Center of the coordination shells.
  //! \param[in] _tolerance : criteria to judge when a shell ends.
  //! \param[in] _natoms : Total number of neighbors to consider. Defaults to fcc + some security.
  PyObject* coordination_shells( crystal::Structure const &_structure, Py_ssize_t _nshells, 
                                 math::rVector3d const &_center,
                                 types::t_real _tolerance=types::tolerance,
                                 Py_ssize_t _natoms = 0 )
    LADA_END({ return (*(PyObject*(*)( crystal::Structure const &, Py_ssize_t,
                                       math::rVector3d const &, types::t_real, 
                                       Py_ssize_t ))
                api_capsule[LADA_SLOT(crystal)])(_structure, _nshells, _center, _tolerance, _natoms); })
#else
  api_capsule[LADA_SLOT(crystal)] = (void *)coordination_shells;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(crystal))
#include LADA_ASSIGN_SLOT(crystal)

#if LADA_CRYSTAL_MODULE != 1
    }
  } // namespace crystal
} // namespace LaDa
#endif
