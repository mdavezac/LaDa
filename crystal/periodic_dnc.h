#if LADA_CRYSTAL_MODULE != 1
namespace LaDa 
{
  namespace crystal 
  {
    namespace
    {
#endif

#if LADA_CRYSTAL_MODULE != 1
  //! Creates divide and conquer box with periodic boundary condition.
  PyObject* dnc_boxes( const Structure &_structure, 
                       math::iVector3d const &_mesh, 
                       types::t_real _overlap)
  LADA_END( { return (*(PyObject*(*)( Structure const&,
                                      math::iVector3d const&,
                                      types::t_real ))
                      api_capsule[BOOST_PP_SLOT(1)])
                     (_structure, _mesh, _overlap); } )
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)dnc_boxes;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
    }
  } // namespace crystal
} // namespace LaDa
#endif
