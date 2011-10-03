#ifndef LADA_NUMPY_TYPES
#define LADA_NUMPY_TYPES

#include <iostream>

#include <boost/mpl/int.hpp>
#include <boost/python/object.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/errors.hpp>
#include <numpy/ndarrayobject.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_common.h>

#include <root_exceptions.h>
#include "exceptions.h"


namespace LaDa
{
  namespace math
  {
    namespace numpy
    {
      //! An mpl integer defining the type.
      template<class T> class type;
      
#     ifdef NUMPY_HAS_LONG_DOUBLE
        //! numpy identifier for long doubles.
        template<> struct type<npy_longdouble> : public boost::mpl::int_<NPY_LONGDOUBLE> 
        {
          //! Original type.
          typedef npy_longdouble np_type;
        };
#     endif
      //! numpy identifier for doubles.
      template<> struct type<npy_double> : public boost::mpl::int_<NPY_DOUBLE> 
      {
        //! Original type.
        typedef npy_double np_type;
      };
      //! numpy identifier for float.
      template<> struct type<npy_float> : public boost::mpl::int_<NPY_FLOAT> 
      {
        //! Original type.
        typedef npy_float np_type;
      };
      //! numpy identifier for long long.
      template<> struct type<npy_longlong> : public boost::mpl::int_<NPY_LONGLONG> 
      {
        //! Original type.
        typedef npy_longlong np_type;
      };
      //! numpy identifier for unsigned long.
      template<> struct type<npy_ulonglong> : public boost::mpl::int_<NPY_ULONGLONG> 
      {
        //! Original type.
        typedef npy_ulonglong np_type;
      };
      //! numpy identifier for long.
      template<> struct type<npy_long> : public boost::mpl::int_<NPY_LONG> 
      {
        //! Original type.
        typedef npy_long np_type;
      };
      //! numpy identifier for unsigned long.
      template<> struct type<npy_ulong> : public boost::mpl::int_<NPY_ULONG> 
      {
        //! Original type.
        typedef npy_ulong np_type;
      };
      //! numpy identifier for int.
      template<> struct type<npy_int> : public boost::mpl::int_<NPY_INT> 
      {
        //! Original type.
        typedef npy_int np_type;
      };
      //! numpy identifier for unsigned int.
      template<> struct type<npy_uint> : public boost::mpl::int_<NPY_UINT> 
      {
        //! Original type.
        typedef npy_uint np_type;
      };
      //! numpy identifier for short.
      template<> struct type<npy_short> : public boost::mpl::int_<NPY_SHORT> 
      {
        //! Original type.
        typedef npy_short np_type;
      };
      //! numpy identifier for unsigned short.
      template<> struct type<npy_ushort> : public boost::mpl::int_<NPY_USHORT> 
      {
        //! Original type.
        typedef npy_ushort np_type;
      };
      //! numpy identifier for byte.
      template<> struct type<npy_byte> : public boost::mpl::int_<NPY_BYTE> 
      {
        //! Original type.
        typedef npy_byte np_type;
      };
      //! numpy identifier for unsigned byte..
      template<> struct type<npy_ubyte> : public boost::mpl::int_<NPY_UBYTE> 
      {
        //! Original type.
        typedef npy_ubyte np_type;
      };
#     ifdef NUMPY_HAS_BOOL
        //! numpy identifier for bool.
        template<> struct type<npy_bool> : public boost::mpl::int_<NPY_BOOL> 
        {
          //! Original type.
          typedef npy_bool np_type;
        };
#     else
        //! numpy identifier for bool.
        template<> struct type<bool> : public boost::mpl::int_<NPY_BOOL> 
        {
          //! Original type.
          typedef npy_bool np_type;
        };
#     endif 
      
      //! Returns true if object is float.
      inline bool is_float(PyObject *_obj_ptr)
      {
        int const nptype = PyArray_ObjectType(_obj_ptr, 0);
        return    nptype ==  NPY_FLOAT
               or nptype ==  NPY_DOUBLE
               or nptype ==  NPY_LONGDOUBLE;
      }
      
      //! Returns true if object is complex.
      inline bool is_complex(PyObject *_obj_ptr)
      {
        int const nptype = PyArray_ObjectType(_obj_ptr, 0);
        return    nptype ==  NPY_CDOUBLE; 
      };
      //! Returns true if object is integer.
      inline bool is_integer(PyObject *_obj_ptr)
      {
        int const nptype = PyArray_ObjectType(_obj_ptr, 0);
        return    nptype ==  NPY_INT
               or nptype ==  NPY_UINT
               or nptype ==  NPY_LONG
               or nptype ==  NPY_ULONG
               or nptype ==  NPY_LONGLONG
               or nptype ==  NPY_ULONGLONG
               or nptype ==  NPY_BYTE
               or nptype ==  NPY_UBYTE
               or nptype ==  NPY_SHORT
               or nptype ==  NPY_USHORT;
      }
      //! Returns true if object is boolean.
      inline bool is_bool(PyObject *_obj_ptr) { return PyArray_ObjectType(_obj_ptr, 0) ==  NPY_BOOL; }
      
      //! Returns true if downcasting from PyObject to T. 
      template<class T> bool is_downcasting(PyObject *_obj_ptr) 
      {
        switch(PyArray_ObjectType(_obj_ptr, 0))
        {
          case NPY_FLOAT     : return sizeof(type<npy_float>::np_type) > sizeof(T);
          case NPY_DOUBLE    : return sizeof(type<npy_double>::np_type) > sizeof(T);
          case NPY_LONGDOUBLE: return sizeof(type<npy_longdouble>::np_type) > sizeof(T);
          case NPY_INT       : return sizeof(type<npy_int>::np_type) > sizeof(T);
          case NPY_UINT      : return sizeof(type<npy_uint>::np_type) > sizeof(T);
          case NPY_LONG      : return sizeof(type<npy_long>::np_type) > sizeof(T);
          case NPY_ULONG     : return sizeof(type<npy_ulong>::np_type) > sizeof(T);
          case NPY_LONGLONG  : return sizeof(type<npy_longlong>::np_type) > sizeof(T);
          case NPY_ULONGLONG : return sizeof(type<npy_ulonglong>::np_type) > sizeof(T);
          case NPY_BYTE      : return sizeof(type<npy_byte>::np_type) > sizeof(T);
          case NPY_UBYTE     : return sizeof(type<npy_ubyte>::np_type) > sizeof(T);
          case NPY_SHORT     : return sizeof(type<npy_short>::np_type) > sizeof(T);
          case NPY_USHORT    : return sizeof(type<npy_ushort>::np_type) > sizeof(T);
          default: break;
        };
        python::PyException<error::ValueError>::throw_error("Unknown numpy array type.");
      }
     
      //! Returns python object from array.
      template<class T>
        boost::python::object array_from_ptr(T *_ptr, size_t n)
        {
          npy_intp dims[1] = {n};
          PyObject *result = PyArray_SimpleNewFromData(1, dims, type<T>::value, _ptr);
          if( result == NULL or PyErr_Occurred() != NULL ) 
          {
            if(PyErr_Occurred() == NULL)
              python::PyException<error::internal>::throw_error("Could not create numpy array");
            boost::python::throw_error_already_set();
            return boost::python::object();
          }
          return boost::python::object(boost::python::handle<>(result));
        }
      //! Returns python object from array, with base.
      template<class T>
        boost::python::object array_from_ptr(T *_ptr, size_t n, PyObject *_base)
        {
          npy_intp dims[1] = {n};
          PyObject *result = PyArray_SimpleNewFromData(1, dims, type<T>::value, _ptr);
          if( result == NULL or PyErr_Occurred() != NULL ) 
          {
            if(PyErr_Occurred() == NULL)
              python::PyException<error::internal>::throw_error("Could not create numpy array");
            boost::python::throw_error_already_set();
            return boost::python::object();
          }
          if(_base != Py_None)
          {
            PyArrayObject* array = (PyArrayObject*)result;
            array->base = _base;
            Py_INCREF(_base);
          }
          return boost::python::object(boost::python::handle<>(result));
        }
      //! Returns python object from array.
      template<class T>
        boost::python::object array_from_ptr(T *_ptr, size_t n0, size_t n1)
        {
          npy_intp dims[2] = {n0, n1};
          PyObject *result = PyArray_SimpleNewFromData(2, dims, type<T>::value, _ptr);
          if( result == NULL or PyErr_Occurred() != NULL ) 
          {
            if(PyErr_Occurred() == NULL)
              python::PyException<error::internal>::throw_error("Could not create numpy array");
            boost::python::throw_error_already_set();
            return boost::python::object();
          }
          return boost::python::object(boost::python::handle<>(result));
        }

      inline boost::python::object _create_array( npy_intp ndims, npy_intp dims[],
                                                  int _type, bool is_fortran = false )
      {
        PyObject *result = PyArray_ZEROS(ndims, dims, _type, is_fortran ? 1: 0);
        if( result == NULL or PyErr_Occurred() != NULL ) 
        {
          if(PyErr_Occurred() == NULL)
            python::PyException<error::internal>::throw_error("Could not create numpy array");
          boost::python::throw_error_already_set();
          return boost::python::object();
        }
        return boost::python::object(boost::python::handle<>(result));
      }

      template<int T_TYPE>
        boost::python::object create_array(npy_intp n0, bool is_fortran = false)
          { npy_intp a[1] = {n0}; return _create_array(1, a, T_TYPE, is_fortran); }
      template<int T_TYPE>
        boost::python::object create_array(npy_intp n0, npy_intp n1, bool is_fortran = false)
          { npy_intp a[2] = {n0, n1}; return _create_array(2, a, T_TYPE, is_fortran); }
      template<int T_TYPE>
        boost::python::object create_array(npy_intp n0, npy_intp n1, npy_intp n2, bool is_fortran = false)
          { npy_intp a[3] = {n0, n1, n2}; return _create_array(3, a, T_TYPE, is_fortran); }
      template<int T_TYPE>
        boost::python::object create_array( npy_intp n0, npy_intp n1, npy_intp n2,
                                                   npy_intp n3, bool is_fortran = false )
          { npy_intp a[4] = {n0, n1, n2, n3}; return _create_array(4, a, T_TYPE, is_fortran); }

     
      inline PyArrayObject* get_pyarray_pointer(boost::python::object const &_object)
         { return reinterpret_cast<PyArrayObject*>(_object.ptr()); }
      template<class T_TYPE>
        inline T_TYPE get_data_pointer(boost::python::object const &_object)
          { return reinterpret_cast<T_TYPE>( reinterpret_cast<PyArrayObject*>(_object.ptr())->data ); }
      inline npy_intp *get_strides_pointer(boost::python::object const &_object)
        { return reinterpret_cast<PyArrayObject*>(_object.ptr())->strides; }

      inline bool check_is_array(boost::python::object const &_ob)
      {
        if( PyArray_Check(_ob.ptr()) ) return true;
        
        python::PyException<error::ValueError>::throw_error("Argument is not a numpy array.\n");
        return false;
      }
      inline bool check_is_complex_array(boost::python::object const &_ob)
      {
        if(not check_is_array(_ob)) return  false;
        if(math::numpy::is_complex(_ob.ptr())) return true;
        
        python::PyException<error::ValueError>::throw_error("Argument is not a complex numpy array.\n");
        return false;
      }

      template<class T> 
        inline boost::python::object copy_1darray(T const &_array)
        {
          int const value = type<typename T::value_type>::value;
          typedef typename type<typename T::value_type>::np_type np_type;
          boost::python::object result = create_array<value>(_array.size());
          PyObject *i_in  = PyArray_IterNew(result.ptr());
          typename T::const_iterator i_first = _array.begin();
          typename T::const_iterator const i_end = _array.end();
          for(; i_first != i_end; ++i_first)
          {
            np_type * const data =  (np_type *)(PyArray_ITER_DATA(i_in));
            *data = *i_first;
            PyArray_ITER_NEXT(i_in);
          }
          Py_DECREF(i_in);
          return result;
        }


      inline boost::python::object from_template( boost::python::object const &_object, 
                                                  int new_type = -1 )
      {
         PyArrayObject *ptr_array = get_pyarray_pointer(_object);
         if(new_type == -1) new_type = PyArray_TYPE(ptr_array);
         bool const fortran = PyArray_ISFORTRAN(ptr_array->ob_type);
         return _create_array(ptr_array->nd, ptr_array->dimensions, new_type, fortran);
      }

      //! Converts data from numpy iterator to requested type.
      template<class T>
        T to_type(PyObject* _iter, int const _type, int const _offset = 0)
        {
          char const * const ptr_data = ((char const *) PyArray_ITER_DATA(_iter) + _offset);
          switch(_type)
          {
            case NPY_FLOAT     :
              return static_cast<T>(*((type<npy_float>::np_type*)      ptr_data));
            case NPY_DOUBLE    :
              return static_cast<T>(*((type<npy_double>::np_type*)     ptr_data));
            case NPY_LONGDOUBLE:
              return static_cast<T>(*((type<npy_longdouble>::np_type*) ptr_data));
            case NPY_INT       :
              return static_cast<T>(*((type<npy_int>::np_type*)        ptr_data));
            case NPY_UINT      :
              return static_cast<T>(*((type<npy_uint>::np_type*)       ptr_data));
            case NPY_LONG      :
              return static_cast<T>(*((type<npy_long>::np_type*)       ptr_data));
            case NPY_ULONG     :
              return static_cast<T>(*((type<npy_ulong>::np_type*)      ptr_data));
            case NPY_LONGLONG  :
              return static_cast<T>(*((type<npy_longlong>::np_type*)   ptr_data));
            case NPY_ULONGLONG :
              return static_cast<T>(*((type<npy_ulonglong>::np_type*)  ptr_data));
            case NPY_BYTE      :
              return static_cast<T>(*((type<npy_byte>::np_type*)       ptr_data));
            case NPY_UBYTE     :
              return static_cast<T>(*((type<npy_ubyte>::np_type*)      ptr_data));
            case NPY_SHORT     :
              return static_cast<T>(*((type<npy_short>::np_type*)      ptr_data));
            case NPY_USHORT    :
              return static_cast<T>(*((type<npy_ushort>::np_type*)     ptr_data));
            default: break;
          };
          python::PyException<error::ValueError>::throw_error("Unknown numpy type on input.");
          return T(-1);
        }
    }
  }
}
#endif

