//
//  Version: $Id$
//
#ifndef _DEBUG_H_
#define _DEBUG_H_

#include <opt/types.h>

//! \brief returns a string with random number
//! \details allows the comparison of two revisions while making sure they are
//! at the same point within the pseudo-random sequence.
std::string keycheck( types::t_int _s=20 );
#endif
