//
//  Version: $Id$
//

#ifndef _ATAT_LAMBDA_IMPL_H_
#  if defined( __VTYPE__ ) || defined (__STYPE__) || defined( __MTYPE__)
#   error Macros __VTYPE__, __STYPE__, or __MTYPE__ already defined.
#  endif
#  define _ATAT_LAMBDA_IMPL_H_ 0
#  define __VTYPE__ LaDa::Eigen::Vector3d
#  define __STYPE__ LaDa::types::t_real
#  define __MTYPE__ LaDa::Eigen::Matrix3d
#elif _ATAT_LAMBDA_IMPL_H_ == 0
#  undef _ATAT_LAMBDA_IMPL_H_
#  undef __VTYPE__ 
#  undef __STYPE__ 
#  undef __MTYPE__ 
#  define _ATAT_LAMBDA_IMPL_H_ 1
#  define __VTYPE__ LaDa::atat::rVector2d
#  define __STYPE__ LaDa::types::t_real
#  define __MTYPE__ LaDa::atat::rMatrix2d
#elif _ATAT_LAMBDA_IMPL_H_ == 1
#  undef _ATAT_LAMBDA_IMPL_H_
#  undef __VTYPE__ 
#  undef __STYPE__ 
#  undef __MTYPE__ 
#  define _ATAT_LAMBDA_IMPL_H_ 2
#  define __VTYPE__ LaDa::Eigen::Vector3i
#  define __STYPE__ LaDa::types::t_int
#  define __MTYPE__ LaDa::Eigen::Matrix3i
#elif _ATAT_LAMBDA_IMPL_H_ == 2
#  undef _ATAT_LAMBDA_IMPL_H_
#  undef __VTYPE__ 
#  undef __STYPE__ 
#  undef __MTYPE__ 
#  define _ATAT_LAMBDA_IMPL_H_ 3
#  define __VTYPE__ LaDa::atat::iVector2d
#  define __STYPE__ LaDa::types::t_int
#  define __MTYPE__ LaDa::atat::iMatrix2d
#else
#  undef _ATAT_LAMBDA_IMPL_H_
#  define _ATAT_LAMBDA_IMPL_H_ 5
#  undef __VTYPE__ 
#  undef __STYPE__ 
#  undef __MTYPE__ 
#endif

#if _ATAT_LAMBDA_IMPL_H_ < 4
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
    struct plain_return_type_2<arithmetic_action<multiply_action>, __VTYPE__, __VTYPE__> {
      typedef __STYPE__ type;
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
