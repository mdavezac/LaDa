#if LADA_CRYSTAL_MODULE != 1
namespace LaDa
{
  namespace crystal
  {
    namespace
    {
#endif

#if LADA_CRYSTAL_MODULE != 1
  //! \brief Creates list of first neigbors up to given input.
  //! \details Always return all nth neighbors. In other words, in fcc, if you
  //!          ask for 6 neighbor, actually 12 are returned. 
  //! \returns A python list of 3-tuples. Each tuple consists of a (unwrapped
  //!          python) reference to an atom, a vector which goes from the
  //!          center to the relevant periodic image of that neighbor, and the
  //!          distance between the center and that neighbor.
  PyObject* neighbors( Structure const &_structure, Py_ssize_t _nmax, 
                       math::rVector3d const &_center,
                       types::t_real _tolerance=types::tolerance )
    LADA_END( { return (*(PyObject*(*)( Structure const&, 
                                        Py_ssize_t, 
                                        math::rVector3d const&, 
                                        types::t_real ))
                        api_capsule[BOOST_PP_SLOT(1)])(_structure, _nmax, _center, _tolerance); } )
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)neighbors;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
    } // anonymous namespace
  } // end of crystal namespace.
} // namespace LaDa
#endif
