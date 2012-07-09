// Some macro to share code between c++ and python functions.
#ifndef LADA_CRYSTAL_HFTRANSFORM_MACRO_HPP
# define LADA_CRYSTAL_HFTRANSFORM_MACRO_HPP
# ifdef LADA_HFTRANSFORM_SHARED0
#   error LADA_HFTRANSFORM_SHARED0 already defined.
# endif
# ifdef LADA_HFTRANSFORM_SHARED1
#   error LADA_HFTRANSFORM_SHARED1 already defined.
# endif
# ifdef LADA_HFTRANSFORM_SHARED2
#   error LADA_HFTRANSFORM_SHARED2 already defined.
# endif
# define LADA_HFTRANSFORM_SHARED0(QUOTIENT, INDEX, SITE)                                      \
    int const flat_result = SITE == -1 ?                                                         \
           INDEX(2) + QUOTIENT(2) * (INDEX(1) + QUOTIENT(1)):                                    \
           INDEX(2) + QUOTIENT(2) * (INDEX(1) + QUOTIENT(1) * (INDEX(0) + SITE * QUOTIENT(0))); 
  
# define LADA_HFTRANSFORM_SHARED1(VECTOR, MATRIX, POS, ERROR, RETURN)   \
       math::iVector3d vector_result;                                      \
       const math::rVector3d pos_(MATRIX*POS);                             \
       const math::iVector3d int_pos                                       \
       (                                                                   \
         types::t_int( rint( pos_(0) ) ),                                  \
         types::t_int( rint( pos_(1) ) ),                                  \
         types::t_int( rint( pos_(2) ) )                                   \
       );                                                                  \
       for( size_t i(0); i < 3; ++i )                                      \
       {                                                                   \
         if( math::neq(pos_(i), types::t_real(int_pos(i))) )               \
         {                                                                 \
           ERROR(ValueError, "Position is not on the lattice.");           \
           RETURN;                                                         \
         }                                                                 \
         vector_result(i) = int_pos(i) % VECTOR(i);                        \
         if( vector_result(i) < 0 ) vector_result(i) += VECTOR(i);         \
       }
# define LADA_HFTRANSFORM_SHARED2(VECTOR) VECTOR(0)*VECTOR(1)*VECTOR(2)
#else 
# undef LADA_HFTRANSFORM_SHARED0
# undef LADA_HFTRANSFORM_SHARED1
# undef LADA_HFTRANSFORM_SHARED2
# undef LADA_CRYSTAL_HFTRANSFORM_MACRO_HPP
#endif
