//
//  Version: $Id$
//
#ifdef _TYPE_
#  undef _TYPE_
#endif
#ifdef _DIM_
#  undef _DIM_
#endif
#ifdef _CLASSNAME_
#  undef _SHORT_TYPE_ 
#endif
#ifdef _PYTHONNAME_
#  undef _PYTHONNAME_ 
#endif
#ifdef _CONCAT_
#  undef _CONCAT_
#endif
#ifdef _CLASSNAME_
#  undef _CLASSNAME_
#endif

#define _CONCAT_( a, b, c, d ) a ## b ## c ## d

#ifndef _WASINCLUDED_
#  define _WASINCLUDED_ 0
#  define _TYPE_ types::t_int 
#  define _DIM_ 2
#  define _CLASSNAME_(object) _CONCAT_( i, object, 2, d )
#  define _PYTHONNAME_(object) "i"#object"2d"
#elif _WASINCLUDED_ == 0
#  undef _WASINCLUDED_ 
#  define _WASINCLUDED_ 1
#  define _TYPE_ types::t_int 
#  define _DIM_  3
#  define _CLASSNAME_(object) _CONCAT_( i, object, 3, d )
#  define _PYTHONNAME_(object) "i"#object"3d"
#elif _WASINCLUDED_ == 1
#  undef _WASINCLUDED_ 
#  define _WASINCLUDED_ 2
#  define _TYPE_ types::t_real 
#  define _DIM_ 2
#  define _CLASSNAME_(object) _CONCAT_( r, object, 2, d )
#  define _PYTHONNAME_(object) "r"#object"2d"
#elif _WASINCLUDED_ == 2
#  undef _WASINCLUDED_ 
#  define _WASINCLUDED_ 3
#  define _TYPE_ types::t_real 
#  define _DIM_  3
#  define _CLASSNAME_(object) _CONCAT_( r, object, 3, d )
#  define _PYTHONNAME_(object) "r"#object"3d"
#else
#  undef _WASINCLUDED_
#  undef __OPP__
#endif

#if defined(_WASINCLUDED_)

#ifndef _INMODULE_

namespace atat
{
  namespace details
  {
#ifndef __OPP__
#define __OPP__
    template < class A, class B, class R > 
    R* mul( const A& _a, const B& _b )
    {
       R* result = new R;
       *result = _a * _b;
       return result;
    }
    template < class A, class B, class R > 
    R* rmul( const A& _a, const B& _b )
      { return mul< B, A, R>( _b, _a); }
    template < class A, class B, class R > 
    R* div( const A& _a, const B& _b )
    {
       R* result = new R;
       *result = _a / _b;
       return result;
    }
    template < class A, class B > 
    A* add( const A& _a, const B& _b )
    {
       A* result = new A;
       *result = _a + _b;
       return result;
    }
    template < class A, class B > 
    B* radd( const A& _a, const B& _b )
      { return add( _b, _a); }
    template < class A, class B > 
    A* diff( const A& _a, const B& _b )
    {
       A* result = new A;
       *result = _a - _b;
       return result;
    }
    template < class A, class B > 
    B* rdiff( const A& _a, const B& _b )
      { return diff( _b, _a); }

    template< types::t_unsigned d > types::t_unsigned length() { return d; }
    template< class T_MAT > T_MAT* transpose( T_MAT &_mat )
      { return new T_MAT( ~_mat ); }
    template< class T_MAT > T_MAT* inverse( T_MAT &_mat )
      { return new T_MAT( !_mat ); }
#endif

    _TYPE_ _CLASSNAME_( getVecItem )( const _CLASSNAME_(Vector) &_v,
                                      types::t_int _i )
    {
      if( _i > _DIM_ or _i < -_DIM_ )
      {
        std::ostringstream sstr;
        sstr << _PYTHONNAME_(Vector) << " holds exactly " 
             << (types::t_unsigned) _DIM_ << " elements. " 
             << "Requested element " << _i << ".";
        throw std::out_of_range( sstr.str() ); 
      }
      return _v.x[ _i < 0 ? _DIM_ + _i: _i ];
    }
    void _CLASSNAME_( setVecItem )( _CLASSNAME_(Vector) &_v,
                                    types::t_int _i, _TYPE_ _a )
    {
      if( _i > _DIM_ or _i < -_DIM_ )
      {
        std::ostringstream sstr;
        sstr << _PYTHONNAME_(Vector) << " holds exactly " 
             << (types::t_unsigned) _DIM_ << " elements. " 
             << "Requested element " << _i << ".";
        throw std::out_of_range( sstr.str() ); 
      }
      _v.x[ _i < 0 ? _DIM_ + _i: _i ] = _a;
    }

    _TYPE_ _CLASSNAME_( getMatItem )( const _CLASSNAME_(Matrix) &_v,
                                      tuple _t )
    {
      types::t_int _i = extract< types::t_int >( _t[0] );
      types::t_int _j = extract< types::t_int >( _t[1] );
      if( _i > _DIM_ or _i < -_DIM_ or
          _j > _DIM_ or _j < -_DIM_ )
      {
        std::ostringstream sstr;
        sstr << _PYTHONNAME_(Matrix) << " holds exactly " 
             << (types::t_unsigned) _DIM_ << "x" 
             << (types::t_unsigned) _DIM_ << " elements. " 
             << "Requested element " << _i << ", " << _j << ".";
        throw std::out_of_range( sstr.str() ); 
      }
      return _v.x[ _i < 0 ? _DIM_ + _i: _i ][ _j < 0 ? _DIM_ + _j: _j ];
    }
    void _CLASSNAME_( setMatItem )( _CLASSNAME_(Matrix) &_v,
                                    tuple _t, _TYPE_ _a )
    {
      types::t_int _i = extract< types::t_int >( _t[0] );
      types::t_int _j = extract< types::t_int >( _t[1] );
      if( _i > _DIM_ or _i < -_DIM_ or
          _j > _DIM_ or _j < -_DIM_ )
      {
        std::ostringstream sstr;
        sstr << _PYTHONNAME_(Vector) << " holds exactly " 
             << (types::t_unsigned) _DIM_ << "x" 
             << (types::t_unsigned) _DIM_ << " elements. " 
             << "Requested element " << _i << ", " << _j << ".";
        throw std::out_of_range( sstr.str() ); 
      }
      _v.x[ _i < 0 ? _DIM_ + _i: _i ][ _j < 0 ? _DIM_ + _j: _j ] = _a;
    }


    std::string _CLASSNAME_(printVector)(const _CLASSNAME_(Vector) &_v)
    { 
      std::ostringstream sstr;
      sstr << _v;
      return sstr.str();
    }
    _CLASSNAME_(Vector)* _CLASSNAME_(makeVector)( object &_ob )
    {
      _CLASSNAME_(Vector)* result = new _CLASSNAME_(Vector);
      for( types::t_int i=_DIM_-1; i >= 0; --i )
        result->x[i] = extract<_TYPE_>(_ob[i]);
      return result;
    }
  }

  _CLASSNAME_(Vector) operator*( const _TYPE_ _b, const _CLASSNAME_(Vector) &_v)
    { return operator*( _v, _b ); }
  bool operator==( const _CLASSNAME_(Vector) &_v, const _CLASSNAME_(Vector) &_b )
    { return std::equal( _v.x, _v.x + _DIM_, _b.x ); }
  bool operator!=( const _CLASSNAME_(Vector) &_v, const _CLASSNAME_(Vector) &_b )
    { return std::equal( _v.x, _v.x + _DIM_, _b.x, 
                         boost::lambda::_1 != boost::lambda::_2 ); }



  bool operator==( const _CLASSNAME_(Matrix) &_v, const _CLASSNAME_(Matrix) &_b )
  {
    for( types::t_int i = 0; i < _DIM_; ++i )
      for( types::t_int j = 0; j < _DIM_; ++j )
        if( _v.x[i][j] != _b.x[i][j] ) return false;
    return true;
  }
  bool operator!=( const _CLASSNAME_(Matrix) &_v, const _CLASSNAME_(Matrix) &_b )
  {
    for( types::t_int i = 0; i < _DIM_; ++i )
      for( types::t_int j = 0; j < _DIM_; ++j )
        if( _v.x[i][j] == _b.x[i][j] ) return false;
    return true;
  }

  _CLASSNAME_(Matrix) operator*( const _TYPE_ _b, const _CLASSNAME_(Matrix) &_v )
    { return operator*(_v, _b); }

  namespace details
  {
    std::string _CLASSNAME_(printMatrix)(const _CLASSNAME_(Matrix) &_v)
    { 
      std::ostringstream sstr;
      sstr << _v;
      return sstr.str();
    }
    _CLASSNAME_(Matrix)* _CLASSNAME_(makeMatrix)( object &_ob )
    {
      _CLASSNAME_(Matrix)* result = new _CLASSNAME_(Matrix);
      for( types::t_int i=_DIM_-1; i >= 0; --i )
        for( types::t_int j=_DIM_-1; j >= 0; --j )
          result->x[i][j] = extract<_TYPE_>(_ob[i][j]);
      return result;
    }
    
  }

}

#else

   class_< FixedVector<_TYPE_, _DIM_> >( _PYTHONNAME_(FixedVector), no_init );
   class_< _CLASSNAME_(Vector), bases< FixedVector<_TYPE_, _DIM_> > >
         ( _PYTHONNAME_(detailsVector), no_init )
     .def( "__add__", 
           &details::add< _CLASSNAME_(Vector),
                          _CLASSNAME_(Vector) >,
           return_value_policy<manage_new_object>() )
     .def( "__diff__", 
           &details::add< _CLASSNAME_(Vector),
                          _CLASSNAME_(Vector) >,
           return_value_policy<manage_new_object>() )
     .def( "__mul__", 
           &details::mul< _CLASSNAME_(Vector),
                          _TYPE_, _CLASSNAME_(Vector) >,
           return_value_policy<manage_new_object>() )
     .def( "__rmul__", 
           &details::rmul< _CLASSNAME_(Vector),
                           _TYPE_, _CLASSNAME_(Vector) >,
           return_value_policy<manage_new_object>() )
     .def( "__div__", 
           &details::div< _CLASSNAME_(Vector),
                          _TYPE_, _CLASSNAME_(Vector) >,
           return_value_policy<manage_new_object>() )
     .def( self == other<_CLASSNAME_(Vector)>() )
     .def( self != other<_CLASSNAME_(Vector)>() )
     .def( "__getitem__", &details::_CLASSNAME_(getVecItem) )
     .def( "__setitem__", &details::_CLASSNAME_(setVecItem) ) 
     .def( "__len__", &details::length<_DIM_> ) 
     .def( "__str__", &details::_CLASSNAME_(printVector) );

   def( _PYTHONNAME_(Vector),
        &details::_CLASSNAME_(makeVector),
        return_value_policy<manage_new_object>() );

   class_< _CLASSNAME_(Matrix) >( _PYTHONNAME_(Matrix) )
     .def( "__mul__", 
           &details::mul< _CLASSNAME_(Matrix),
                          _CLASSNAME_(Matrix),
                          _CLASSNAME_(Matrix) >,
           return_value_policy<manage_new_object>() )
     .def( "__add__", 
           &details::add< _CLASSNAME_(Matrix),
                          _CLASSNAME_(Matrix) >,
           return_value_policy<manage_new_object>() )
     .def( "__diff__", 
           &details::diff< _CLASSNAME_(Matrix),
                          _CLASSNAME_(Matrix) >,
           return_value_policy<manage_new_object>() )
     .def( "__mul__", 
           &details::mul< _CLASSNAME_(Matrix),
                          _TYPE_, _CLASSNAME_(Matrix) >,
           return_value_policy<manage_new_object>() )
     .def( "__rmul__", 
           &details::rmul< _CLASSNAME_(Matrix),
                           _TYPE_, _CLASSNAME_(Matrix) >,
           return_value_policy<manage_new_object>() )
     .def( "__div__", 
           &details::mul< _CLASSNAME_(Matrix),
                          _TYPE_, _CLASSNAME_(Matrix) >,
           return_value_policy<manage_new_object>() )
     .def( "__mul__", 
           &details::mul< _CLASSNAME_(Matrix),
                          _CLASSNAME_(Vector),
                          _CLASSNAME_(Vector) >,
           return_value_policy<manage_new_object>() )
     .def( self == other<_CLASSNAME_(Matrix)>() )
     .def( self != other<_CLASSNAME_(Matrix)>() )
     .def( "__str__", &details::_CLASSNAME_(printMatrix) )
     .def( "__getitem__", &details::_CLASSNAME_(getMatItem) )
     .def( "__setitem__", &details::_CLASSNAME_(setMatItem) ) 
     .def( "__len__", &details::length<_DIM_> ) 
     .def( "transpose", details::transpose< _CLASSNAME_(Matrix) >,
           return_value_policy< manage_new_object >() )
     .def( "inverse", details::inverse< _CLASSNAME_(Matrix) >,
           return_value_policy< manage_new_object >() )
     .def( "diag", &_CLASSNAME_(Matrix)::diag );

   def( "make_"_PYTHONNAME_(Matrix),
        &details::_CLASSNAME_(makeMatrix),
        return_value_policy<manage_new_object>() );

#endif

#include "atat.py.hpp"

#endif 
