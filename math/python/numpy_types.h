#ifndef LADA_NUMPY_TYPES
#define LADA_NUMPY_TYPES

#include <boost/mpl/int.hpp>
#include <numpy/ndarrayobject.h>
#include <numpy/npy_common.h>

namespace LaDa
{
  namespace math
  {
    namespace numpy
    {
      //! An mpl integer defining the type.
      template<class T> class type;
      
      //! numpy identifier for long doubles.
      template<> struct type<npy_longdouble> : public boost::mpl::int_<NPY_LONGDOUBLE> 
      {
        //! Original type.
        typedef npy_longdouble np_type;
      };
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
        typedef npu_longlong np_type;
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
      //! numpy identifier for bool.
      template<> struct type<npy_bool> : public boost::mpl::int_<NPY_BOOL> 
      {
        //! Original type.
        typedef npy_ubyte np_type;
      };
      
      //! Returns true if object is float.
      bool is_float(PyObject *_obj_ptr)
      {
        int const nptype = PyArray_ObjectType(_obj_ptr, 0);
        return    npy_type ==  NPY_FLOAT
               or npy_type ==  NPY_DOUBLE
               or npy_type ==  NPY_LONGDOUBLE;
      }
      
      //! Returns true if object is integer.
      bool is_integer(PyObject *_obj_ptr)
      {
        int const nptype = PyArray_ObjectType(_obj_ptr, 0);
        return    npy_type ==  NPY_INT
               or npy_type ==  NPY_UINT
               or npy_type ==  NPY_LONG
               or npy_type ==  NPY_ULONG
               or npy_type ==  NPY_LONGLONG
               or npy_type ==  NPY_ULONGLONG
               or npy_type ==  NPY_BYTE
               or npy_type ==  NPY_UBYTE
               or npy_type ==  NPY_SHORT
               or npy_type ==  NPY_USHORT
      }
      //! Returns true if object is boolean.
      bool is_integer(PyObject *_obj_ptr) { return PyArray_ObjectType(_obj_ptr, 0) ==  NPY_BOOL; }
      
      //! Returns true if downcasting from PyObject to T. 
      template<class T> bool is_downcasting(PyObject *_obj_ptr) 
      {
        int const nptype = PyArray_ObjectType(_obj_ptr, 0);
        switch(nptype)
        {
          case NPY_FLOAT     : return sizeof(type<NPY_FLOAT>::np_type) > sizeof(T);
          case NPY_DOUBLE    : return sizeof(type<NPY_DOUBLE>::np_type) > sizeof(T);
          case NPY_LONGDOUBLE: return sizeof(type<NPY_LONGDOUBLE>::np_type) > sizeof(T);
          case NPY_INT       : return sizeof(type<NPY_INT>::np_type) > sizeof(T);
          case NPY_UINT      : return sizeof(type<NPY_UINT>::np_type) > sizeof(T);
          case NPY_LONG      : return sizeof(type<NPY_LONG>::np_type) > sizeof(T);
          case NPY_ULONG     : return sizeof(type<NPY_ULONG>::np_type) > sizeof(T);
          case NPY_LONGLONG  : return sizeof(type<NPY_LONGLONG>::np_type) > sizeof(T);
          case NPY_ULONGLONG : return sizeof(type<NPY_ULONGLONG>::np_type) > sizeof(T);
          case NPY_BYTE      : return sizeof(type<NPY_BYTE>::np_type) > sizeof(T);
          case NPY_UBYTE     : return sizeof(type<NPY_UBYTE>::np_type) > sizeof(T);
          case NPY_SHORT     : return sizeof(type<NPY_SHORT>::np_type) > sizeof(T);
          case NPY_USHORT    : return sizeof(type<NPY_USHORT>::np_type) > sizeof(T);
          default: break;
        };
        throw;
      }

    }
  }
}
#endif

