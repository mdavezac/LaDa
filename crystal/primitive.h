#if LADA_CRYSTAL_MODULE != 1

namespace LaDa
{
  namespace crystal 
  {
    namespace
    {
#endif

#if LADA_CRYSTAL_MODULE != 1
  //! Returns the primitive unit structure. 
  Structure primitive(Structure const &_structure, types::t_real _tolerance = -1e0)
      LADA_END({ return (*(Structure(*)(Structure const&, types::t_real))
                         api_capsule[BOOST_PP_SLOT(1)])(_structure, _tolerance); })
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)primitive;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
  //! Returns True if the input is primitive.
  bool is_primitive(Structure const &_structure, types::t_real _tolerance = -1e0)
    LADA_END( { return (*(bool(*)(Structure const&, types::t_real))
                        api_capsule[BOOST_PP_SLOT(1)])(_structure, _tolerance); } )
#else
  api_capsule[BOOST_PP_SLOT(1)] = (void *)is_primitive;
#endif
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)

#if LADA_CRYSTAL_MODULE != 1
    } // anonymous namespace
  } // namespace crystal
} // namespace LaDa
#endif
