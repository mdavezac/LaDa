#ifndef LADA_CRYSTAL_PYTHON_ATOM_HPP
#define LADA_CRYSTAL_PYTHON_ATOM_HPP

#include "LaDaConfig.h"

#include <boost/python/stl_iterator.hpp>
#include <boost/python/object.hpp>

#include <math/python/python.hpp>

#include "../../traits.h"
#include "../../is_container.h"

namespace LaDa
{
  namespace python
  {
    //! checks whether this is a specie.
    template<class T> typename boost::enable_if< crystal::details::is_scalar<T>, bool > :: type
      is_specie(boost::python::object const &_in)
        { return boost::python::extract<T>(_in).check(); }
    //! checks whether this is a string specie.
    template<> bool is_specie<std::string>(boost::python::object const &_in)
        { return PyString_Check(_in.ptr()); }
    //! checks whether this is a specie index.
    template<class T> typename boost::enable_if< crystal::details::is_iterable<T>, bool > :: type
      is_specie(boost::python::object const &_in)
      {
        if(boost::python::extract<T&>(_in).check()) return true;
        if(boost::python::extract<typename T::value_type&>(_in).check()) return true;
        if(not PySequence_Check(_in.ptr()))
          PyException<error::TypeError>::throw_error("Input is not a sequence.");
        boost::python::object i_specie( boost::python::handle<>(PyObject_GetIter(_in.ptr())) );
        if(i_specie.ptr() == Py_None) return true;
        while(PyObject *specie = PyIter_Next(i_specie.ptr()))
        {
          if(not boost::python::extract<typename T::value_type>(specie).check())
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
      extract_specie(boost::python::object _in, T &_type)
        { _type = boost::python::extract<T>(_in); }
    //! Extracts set type to c++.
    template<class T> typename boost::enable_if< crystal::details::is_set<T>, void >::type
      extract_specie(boost::python::object _in, T &_type)
      {
        if(boost::python::extract<T&>(_in).check())
          _type = boost::python::extract<T&>(_in);
        else if(boost::python::extract<typename T::value_type&>(_in).check())
          _type.insert( boost::python::extract<typename T::value_type const &>(_in)() );
        else
        {
          boost::python::stl_input_iterator<typename T::value_type> i_first(_in);
          boost::python::stl_input_iterator<typename T::value_type> const i_end;
          for(; i_first != i_end; ++i_first) _type.insert(*i_first);
        }
      }
    //! Extracts vector type to c++.
    template<class T> typename boost::enable_if< crystal::details::is_container<T>, void >::type
      extract_specie(boost::python::object _in, T &_type)
      {
        if(boost::python::extract<T&>(_in).check())
          _type = boost::python::extract<T&>(_in);
        else if(boost::python::extract<typename T::value_type&>(_in).check())
          _type.push_back( boost::python::extract<typename T::value_type&>(_in)() );
        else
        {
          boost::python::stl_input_iterator<typename T::value_type> i_first(_in);
          boost::python::stl_input_iterator<typename T::value_type> const i_end;
          for(; i_first != i_end; ++i_first) _type.push_back(*i_first);
        }
      }
  }
} 
#endif
