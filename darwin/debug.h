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



//! prints a container as though it were a bistring.
template< class T_CONTAINER >
std::string bitprint( const T_CONTAINER &_container )
{
  std::string result; result.resize( _container.size() );
  typename T_CONTAINER :: const_iterator i_first = _container.begin();
  typename T_CONTAINER :: const_iterator i_last = _container.end();
  for(types::t_unsigned i=0; i_first != i_last; ++i_first, ++i)
    result[i] = ( *i_first > 0 ? '1': '0' );
  return result;
}
#endif
