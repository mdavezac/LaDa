#ifndef LADA_MATH_LAMBDA_IMPL_H
#  if defined( __VTYPE__ ) || defined (__STYPE__) || defined( __MTYPE__)
#   error Macros __VTYPE__, __STYPE__, or __MTYPE__ already defined.
#  endif
#  define LADA_MATH_LAMBDA_IMPL_H 0
#  define __VTYPE__ LaDa::math::rVector3d
#  define __STYPE__ LaDa::types::t_real
#  define __MTYPE__ LaDa::math::rMatrix3d
#elif LADA_MATH_LAMBDA_IMPL_H == 0
#  undef LADA_MATH_LAMBDA_IMPL_H
#  undef __VTYPE__ 
#  undef __STYPE__ 
#  undef __MTYPE__ 
#  define LADA_MATH_LAMBDA_IMPL_H 1
#  define __VTYPE__ LaDa::math::iVector3d
#  define __STYPE__ LaDa::types::t_int
#  define __MTYPE__ LaDa::math::iMatrix3d
#elif LADA_MATH_LAMBDA_IMPL_H == 1
#  undef LADA_MATH_LAMBDA_IMPL_H
#  define LADA_MATH_LAMBDA_IMPL_H 2
#  undef __VTYPE__ 
#  undef __STYPE__ 
#  undef __MTYPE__ 
#endif

#if LADA_MATH_LAMBDA_IMPL_H < 2
#include <boost/lambda/control_structures.hpp>
//! \cond
namespace boost { 
  namespace lambda {
    
    template<> 
    struct plain_return_type_2<arithmetic_action<plus_action>, __VTYPE__, __VTYPE__> {
      typedef __VTYPE__ type;
    };
    template<> 
    struct plain_return_type_2<arithmetic_action<minus_action>, __VTYPE__, __VTYPE__> {
      typedef __VTYPE__ type;
    };

    template<> 
    struct plain_return_type_2<arithmetic_action<multiply_action>, __VTYPE__, __STYPE__> {
      typedef __VTYPE__ type;
    };
    template<> 
    struct plain_return_type_2<arithmetic_action<multiply_action>, __STYPE__, __VTYPE__> {
      typedef __VTYPE__ type;
    };
    template<> 
    struct plain_return_type_2<arithmetic_action<divide_action>, __VTYPE__, __STYPE__> {
      typedef __VTYPE__ type;
    };
    template<> 
    struct plain_return_type_2<arithmetic_action<divide_action>, __STYPE__, __VTYPE__> {
      typedef __VTYPE__ type;
    };


    template<> 
    struct plain_return_type_2<arithmetic_action<plus_action>, __MTYPE__, __MTYPE__> {
      typedef __MTYPE__ type;
    };
    template<> 
    struct plain_return_type_2<arithmetic_action<minus_action>, __MTYPE__, __MTYPE__> {
      typedef __MTYPE__ type;
    };
    template<> 
    struct plain_return_type_2<arithmetic_action<multiply_action>, __MTYPE__, __MTYPE__> {
      typedef __MTYPE__ type;
    };
    template<> 
    struct plain_return_type_1<logical_action<not_action>, __MTYPE__> {
      typedef __MTYPE__ type;
    };

    template<> 
    struct plain_return_type_2<arithmetic_action<multiply_action>, __MTYPE__, __VTYPE__> {
      typedef __VTYPE__ type;
    };
    template<> 
    struct plain_return_type_2<arithmetic_action<multiply_action>, __MTYPE__, __STYPE__> {
      typedef __MTYPE__ type;
    };
    template<> 
    struct plain_return_type_2<arithmetic_action<multiply_action>, __STYPE__, __MTYPE__> {
      typedef __MTYPE__ type;
    };
    template<> 
    struct plain_return_type_2<arithmetic_action<divide_action>, __MTYPE__, __STYPE__> {
      typedef __MTYPE__ type;
    };
    template<> 
    struct plain_return_type_2<arithmetic_action<divide_action>, __STYPE__, __MTYPE__> {
      typedef __MTYPE__ type;
    };

  }
}
//! \endcond
// reincludes itself for other types.
#include "lambda.impl.h" 

#endif
