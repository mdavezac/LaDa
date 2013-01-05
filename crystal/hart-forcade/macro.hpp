// Some macro to share code between c++ and python functions.
#ifndef PYLADA_CRYSTAL_HFTRANSFORM_MACRO_HPP
# define PYLADA_CRYSTAL_HFTRANSFORM_MACRO_HPP
# ifdef PYLADA_HFTRANSFORM_SHARED0
#   error PYLADA_HFTRANSFORM_SHARED0 already defined.
# endif
# ifdef PYLADA_HFTRANSFORM_SHARED1
#   error PYLADA_HFTRANSFORM_SHARED1 already defined.
# endif
# ifdef PYLADA_HFTRANSFORM_SHARED2
#   error PYLADA_HFTRANSFORM_SHARED2 already defined.
# endif
# define PYLADA_HFTRANSFORM_SHARED0(QUOTIENT, INDEX, SITE)                                         \
    int const flat_result = SITE == -1 ?                                                         \
           INDEX(2) + QUOTIENT(2) * (INDEX(1) + QUOTIENT(1) * INDEX(0)):                         \
           INDEX(2) + QUOTIENT(2) * (INDEX(1) + QUOTIENT(1) * (INDEX(0) + SITE * QUOTIENT(0))); 
  
# define PYLADA_HFTRANSFORM_SHARED1(VECTOR, MATRIX, POS, ERROR, RETURN)   \
       math::iVector3d vector_result;                                      \
       const math::rVector3d pos_(MATRIX*POS);                             \
       const math::iVector3d int_pos                                       \
       (                                                                   \
         types::t_int( rint( pos_(0) + 1e-8 ) ),                           \
         types::t_int( rint( pos_(1) + 1e-8 ) ),                           \
         types::t_int( rint( pos_(2) + 1e-8 ) )                            \
       );                                                                  \
       for( size_t i(0); i < 3; ++i )                                      \
       {                                                                   \
         if( math::neq(pos_(i), types::t_real(int_pos(i)), 1e-4) )         \
         {                                                                 \
           ERROR(ValueError, "Position is not on the lattice.");           \
           RETURN;                                                         \
         }                                                                 \
         vector_result(i) = int_pos(i) % VECTOR(i);                        \
         if( vector_result(i) < 0 ) vector_result(i) += VECTOR(i);         \
       }
# define PYLADA_HFTRANSFORM_SHARED2(VECTOR) VECTOR(0)*VECTOR(1)*VECTOR(2)
#else 
# undef PYLADA_HFTRANSFORM_SHARED0
# undef PYLADA_HFTRANSFORM_SHARED1
# undef PYLADA_HFTRANSFORM_SHARED2
# undef PYLADA_CRYSTAL_HFTRANSFORM_MACRO_HPP
#endif
