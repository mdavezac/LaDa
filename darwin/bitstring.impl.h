//
//  Version: $Id: bitstring.h 332 2007-10-20 00:44:05Z davezac $
//
#ifndef _BITSTRING_IMPL_H_
#define _BITSTRING_IMPL_H_

namespace BitString
{
  template<class T_CONT>
  std::ostream& operator<<(std::ostream &_stream, const Object<T_CONT> &_o)
  {
    typedef typename Object<T_CONT> ::t_Type t_Type;
    typedef typename Object<T_CONT> ::t_Container :: const_iterator const_iterator;
    const_iterator i_var = _o.bitstring.begin();
    const_iterator i_end = _o.bitstring.end();
    for(; i_var != i_end; ++i_var )
      _stream << ( spin_up<t_Type>(*i_var) ? '1' : '0' );
    return _stream;
  }
  template<class T_CONT>
  void operator<<(std::string &_str, const Object<T_CONT> &_o)
  {
    std::ostringstream sstr;
    sstr << _o; _str = sstr.str();
  }
  template<class T_CONT>
  void operator<<(Object<T_CONT> &_o, const std::string &_c)
  {
    typedef typename Object<T_CONT> ::t_Type t_Type;
    typedef typename Object<T_CONT> ::t_Container :: iterator iterator;

    types::t_unsigned size = _c.size();
    _o.bitstring.resize( size );
    iterator i_var = _o.bitstring.begin();
    iterator i_end = _o.bitstring.end();
    for(types::t_unsigned n=0; i_var != i_end; ++i_var, ++n )
      *i_var = ( spin_up<t_Type>(_c[n]) ? 1.0: -1.0 );
  }


} // namespace BitString

#endif
