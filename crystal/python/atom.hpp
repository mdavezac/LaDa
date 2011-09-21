#ifndef LADA_CRYSTAL_PYTHON_ATOM_HPP
#define LADA_CRYSTAL_PYTHON_ATOM_HPP

#include "LaDaConfig.h"

#include "../traits.h"
#include "../is_container.h"

namespace LaDa
{
  namespace python
  {
    void expose_atom();
    namespace details
    {
      template<class T> typename boost::enable_if< crystal::details::is_scalar<T>, void >::type
        extract_type0(bp::object _in, T &_type)
          { _type = bp::extract<T>(_in); }
      template<class T> typename boost::enable_if< crystal::details::is_scalar<T> , void>::type
        extract_type1(bp::object _in, T &_type)
        {
          if(bp::len(_in) > 1) 
            PyException<error::TypeError>::throw_error("Atomic type needs be a scalar.");
          _type = bp::extract<T>(_in); 
        }
      template<class T> typename boost::enable_if< crystal::details::is_set<T>, void >::type
        extract_type0(bp::object _in, T &_type)
        {
          if(PyObject_HasAttrString(_in.ptr(), "__iter__") == 0)
            _type.insert( bp::extract<typename T::value_type>(_in) );
          else
          {
            bp::stl_input_iterator<typename T::value_type> i_first(_in);
            bp::stl_input_iterator<typename T::value_type> const i_end;
            for(; i_first != i_end; ++i_first) _type.insert(*i_first);
          }
        }
      template<class T> typename boost::enable_if< crystal::details::is_container<T>, void >::type
        extract_type0(bp::object _in, T &_type)
        {
          if(PyObject_HasAttrString(_in.ptr(), "__iter__") == 0)
            _type.push_back( bp::extract<typename T::value_type>(_in) );
          else
          {
            bp::stl_input_iterator<typename T::value_type> i_first(_in);
            bp::stl_input_iterator<typename T::value_type> const i_end;
            for(; i_first != i_end; ++i_first) _type.push_back(*i_first);
          }
        }
      template<class T> typename boost::enable_if< crystal::details::is_iterable<T> , void>::type
        extract_type1(bp::object _in, T &_type)
          { extract_type0(_in, _type); }
    }
  }
} 
#endif
