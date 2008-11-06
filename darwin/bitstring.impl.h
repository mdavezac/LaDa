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
#     ifdef _ALLOY_LAYERS_
        _stream << *i_var << " ";
#     else
        _stream << ( spin_up<t_Type>(*i_var) ? '1' : '0' );
#     endif
    return _stream;
  }
  template<class T_CONT>
  void operator<<(std::string &_str, const Object<T_CONT> &_o)
  {
    std::ostringstream sstr;
    typedef typename Object<T_CONT> ::t_Type t_Type;
    typedef typename Object<T_CONT> ::t_Container :: const_iterator const_iterator;
    const_iterator i_var = _o.bitstring.begin();
    const_iterator i_end = _o.bitstring.end();
    for(; i_var != i_end; ++i_var )
#     ifdef _ALLOY_LAYERS_
        sstr << *i_var << " ";
#     else
        sstr << ( spin_up<t_Type>(*i_var) ? '1' : '0' );
#     endif
    _str = sstr.str();
  }
  template<class T_CONT>
  void operator<<(Object<T_CONT> &_o, const std::string &_c)
  {
#   ifdef _ALLOY_LAYERS_
    std::istringstream sstr( _c );
    _o.bitstring.clear();
    while( sstr.good() )
    {
      typedef typename T_CONT :: value_type t_Type;
      t_Type bit;
      sstr >> bit;
      _o.bitstring.push_back( bit );
    }
#   else
      typedef typename Object<T_CONT> ::t_Type t_Type;
      typedef typename Object<T_CONT> ::t_Container :: iterator iterator;
  
      types::t_unsigned size = _c.size();
      _o.bitstring.resize( size );
      iterator i_var = _o.bitstring.begin();
      iterator i_end = _o.bitstring.end();
      for(types::t_unsigned n=0; i_var != i_end; ++i_var, ++n )
        *i_var = ( spin_up<t_Type>(_c[n]) ? 1.0: -1.0 );
#   endif
  }


} // namespace BitString

#endif
