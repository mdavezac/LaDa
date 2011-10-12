#ifndef LADA_CRYSTAL_PYTHON_ATOM_HPP
#define LADA_CRYSTAL_PYTHON_ATOM_HPP

#include "LaDaConfig.h"

#include <math/python/python.hpp>

#include "../traits.h"
#include "../is_container.h"

namespace LaDa
{
  namespace python
  {
    //! checks whether this is a specie.
    template<class T> typename boost::enable_if< crystal::details::is_scalar<T>, bool > :: type
      is_specie(bp::object const &_index)
        { return bp::extract<T>(_index).check(); }
    //! checks whether this is a string specie.
    template<> typename bool is_specie(bp::object const &_index)
        { return PyString_Check(_index.ptr()); }
    //! checks whether this is a specie index.
    template<class T> typename boost::enable_if< crystal::details::is_iterable<T>, bool > :: type
      is_specie(bp::object const &_index)
      {
        if(bp::extract<T&>(_index).check()) return true;
        if(bp::extract<typename T::value_type&>(_index).check()) return true;
        if(not PySequence_Check(_index.ptr()))
          PyException<error::TypeError>::throw_error("Input is not a sequence.");
        bp::object i_specie( bp::handle<>(PyObject_GetIter(_index.ptr())) );
        if(i_specie.ptr() == Py_None) return true;
        while(PyObject *specie = PyIter_Next(i_specie.ptr()))
        {
          if(not bp::extract<typename T::value_type>(specie).check())
          {
            Py_DECREF(specie);
            return false;
          }
          Py_DECREF(specie);
        }
        return true;
      };
    //! Extracts scalar type to c++.
    template<class T> typename boost::enable_if< crystal::details::is_scalar<T>, void >::type
      extract_specie(bp::object _in, T &_type)
        { _type = bp::extract<T>(_in); }
    //! Extracts set type to c++.
    template<class T> typename boost::enable_if< crystal::details::is_set<T>, void >::type
      extract_specie(bp::object _in, T &_type)
      {
        if(bp::extract<T&>(_in).check())
          _type = bp::extract<T&>(_in);
        else if(bp::extract<typename T::value_type&>(_index).check())
          _type.insert( bp::extract<typename T::value_type const &>(_in)() );
        else
        {
          bp::stl_input_iterator<typename T::value_type> i_first(_in);
          bp::stl_input_iterator<typename T::value_type> const i_end;
          for(; i_first != i_end; ++i_first) _type.insert(*i_first);
        }
      }
    //! Extracts vector type to c++.
    template<class T> typename boost::enable_if< crystal::details::is_container<T>, void >::type
      extract_specie(bp::object _in, T &_type)
      {
        if(bp::extract<T&>(_in).check())
          _type = bp::extract<T&>(_in);
        else if(bp::extract<typename T::value_type&>(_index).check())
          _type.push_back( bp::extract<typename T::value_type&>(_in)() );
        else
        {
          bp::stl_input_iterator<typename T::value_type> i_first(_in);
          bp::stl_input_iterator<typename T::value_type> const i_end;
          for(; i_first != i_end; ++i_first) _type.push_back(*i_first);
        }
      }
  }
} 
#endif
