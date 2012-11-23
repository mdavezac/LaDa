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
                      api_capsule[LADA_SLOT(crystal)])
                     (_structure, _mesh, _overlap); } )
#else
  api_capsule[LADA_SLOT(crystal)] = (void *)dnc_boxes;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(LADA_SLOT(crystal))
#include LADA_ASSIGN_SLOT(crystal)

#if LADA_CRYSTAL_MODULE != 1
    }
  } // namespace crystal
} // namespace LaDa
#endif
