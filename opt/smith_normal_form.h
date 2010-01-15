//
//  Version: $Id$
//

#ifndef _LADA_OPT_SMITH_NORMAL_FORM_H_
#define _LADA_OPT_SMITH_NORMAL_FORM_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Eigen/Core> 

namespace LaDa 
{
  namespace opt
  {
    //! Computes smith normal form of a matrix \a  _S = \a _L \a _M \a _R.
    void smith_normal_form( Eigen::Matrix3i& _S, Eigen::Matrix3i & _L,
                            const Eigen::Matrix3i& _M, Eigen::Matrix3i &_R );
  } // namespace opt

} // namespace LaDa

#endif
