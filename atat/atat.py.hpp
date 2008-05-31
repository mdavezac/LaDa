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
# undef _WASINCLUDED_
#endif

#if defined(_WASINCLUDED_)

#ifndef _INMODULE_

namespace atat
{

  struct _CLASSNAME_(WVector) : public _CLASSNAME_(Vector)
  {
    typedef _CLASSNAME_(Vector) t_Base;
    public:
      _CLASSNAME_(WVector)() { std::fill( x, x + _DIM_, _TYPE_(0) ); }
      _CLASSNAME_(WVector) ( const _CLASSNAME_(WVector) &_V )
        { std::copy( _V.x, _V.x + _DIM_, x ); }
      _CLASSNAME_(WVector) ( const _CLASSNAME_(Vector) &_V )
        { std::copy( _V.x, _V.x + _DIM_, x ); }
      _CLASSNAME_(WVector) ( const boost::python::object & _ob );
      std::string print() const;
  };

  _CLASSNAME_(WVector)* operator+( const _CLASSNAME_(WVector) &_v, const _CLASSNAME_(WVector) &_b )
  {
    _CLASSNAME_(WVector) *result = new _CLASSNAME_(WVector);
    std::transform( _v.x, _v.x + _DIM_, _b.x, result->x,
                    boost::lambda::_1 + boost::lambda::_2 );
    return result;
  }
  _CLASSNAME_(WVector)* operator-( const _CLASSNAME_(WVector) &_v, const _CLASSNAME_(WVector) &_b )
  {
    _CLASSNAME_(WVector) *result = new _CLASSNAME_(WVector);
    std::transform( _v.x, _v.x + _DIM_, _b.x, result->x,
                    boost::lambda::_1 - boost::lambda::_2 );
    return result;
  }
  _CLASSNAME_(WVector)* operator*( const _CLASSNAME_(WVector) &_v, const _TYPE_ _b )
  {
    _CLASSNAME_(WVector) *result = new _CLASSNAME_(WVector);
    std::transform( _v.x, _v.x + _DIM_, result->x,
                    boost::lambda::_1 * boost::lambda::constant(_b) );
    return result;
  }
  _CLASSNAME_(WVector)* operator/( const _CLASSNAME_(WVector) &_v, const _TYPE_ _b )
  {
    _CLASSNAME_(WVector) *result = new _CLASSNAME_(WVector);
    std::transform( _v.x, _v.x + _DIM_, result->x,
                    boost::lambda::_1 / boost::lambda::constant(_b) );
    return result;
  }
  _CLASSNAME_(WVector)* operator*( const _TYPE_ _b, const _CLASSNAME_(WVector) &_v)
    { return operator*( _v, _b ); }
  bool operator==( const _CLASSNAME_(WVector) &_v, const _CLASSNAME_(WVector) &_b )
    { return std::equal( _v.x, _v.x + _DIM_, _b.x ); }
  bool operator!=( const _CLASSNAME_(WVector) &_v, const _CLASSNAME_(WVector) &_b )
    { return std::equal( _v.x, _v.x + _DIM_, _b.x, 
                         boost::lambda::_1 != boost::lambda::_2 ); }



  _CLASSNAME_(WVector) :: _CLASSNAME_(WVector)( const boost::python::object &_ob )
  {
    for( types::t_int i = 0; i < _DIM_; ++i )
      x[i] = boost::python::extract<_TYPE_>( _ob[i] ); 
  }

  std::ostream& operator<<( std::ostream &_o, const _CLASSNAME_(WVector)& _v )
  {
    std::for_each( _v.x, _v.x + _DIM_, 
                   boost::lambda::var( _o ) << boost::lambda::_1
                                            << boost::lambda::constant(" ") );
    return _o;
  }
  std::string _CLASSNAME_(WVector) :: print() const
  {
    std::ostringstream sstr;
    sstr << *this;
    return sstr.str();
  }



  struct _CLASSNAME_(WMatrix) : public FixedMatrix<_TYPE_, _DIM_> 
  {
    public:
      _CLASSNAME_(WMatrix)() { zero(); }
      _CLASSNAME_(WMatrix) ( const _CLASSNAME_(WMatrix) &_V ) : FixedMatrix<_TYPE_, _DIM_>( _V ) {}
      _CLASSNAME_(WMatrix) ( const boost::python::object & _ob );
      std::string print() const;
      void operator+=( const _CLASSNAME_(WMatrix) &_b );
      void operator-=( const _CLASSNAME_(WMatrix) &_b );
      void diag( const _CLASSNAME_(WVector) & _v );
  };

  void _CLASSNAME_(WMatrix) :: diag( const _CLASSNAME_(WVector) &_v )
  {
    for( types::t_int i = 0; i < _DIM_; ++i )
      for( types::t_int j = 0; j < _DIM_; ++j )
        x[i][j] = i == j ? _v[i]: _TYPE_(0);
  }


  void _CLASSNAME_(WMatrix) :: operator+=( const _CLASSNAME_(WMatrix) &_v )
  {
    for( types::t_int i = 0; i < _DIM_; ++i )
      for( types::t_int j = 0; j < _DIM_; ++j )
        x[i][j] +=  _v.x[i][j];
  }
  void _CLASSNAME_(WMatrix) :: operator-=( const _CLASSNAME_(WMatrix) &_v )
  {
    for( types::t_int i = 0; i < _DIM_; ++i )
      for( types::t_int j = 0; j < _DIM_; ++j )
        x[i][j] +=  _v.x[i][j];
  }
  _CLASSNAME_(WMatrix)* operator+( const _CLASSNAME_(WMatrix) &_v, const _CLASSNAME_(WMatrix) &_b )
  {
    _CLASSNAME_(WMatrix) *result = new _CLASSNAME_(WMatrix);
    for( types::t_int i = 0; i < _DIM_; ++i )
      for( types::t_int j = 0; j < _DIM_; ++j )
        result->x[i][j] =  _v.x[i][j] + _b.x[i][j]; 
    return result;
  }
  _CLASSNAME_(WMatrix)* operator-( const _CLASSNAME_(WMatrix) &_v, const _CLASSNAME_(WMatrix) &_b )
  {
    _CLASSNAME_(WMatrix) *result = new _CLASSNAME_(WMatrix);
    for( types::t_int i = 0; i < _DIM_; ++i )
      for( types::t_int j = 0; j < _DIM_; ++j )
        result->x[i][j] =  _v.x[i][j] - _b.x[i][j]; 
    return result;
  }
  bool operator==( const _CLASSNAME_(WMatrix) &_v, const _CLASSNAME_(WMatrix) &_b )
  {
    for( types::t_int i = 0; i < _DIM_; ++i )
      for( types::t_int j = 0; j < _DIM_; ++j )
        if( _v.x[i][j] != _b.x[i][j] ) return false;
    return true;
  }
  bool operator!=( const _CLASSNAME_(WMatrix) &_v, const _CLASSNAME_(WMatrix) &_b )
  {
    for( types::t_int i = 0; i < _DIM_; ++i )
      for( types::t_int j = 0; j < _DIM_; ++j )
        if( _v.x[i][j] == _b.x[i][j] ) return false;
    return true;
  }

  _CLASSNAME_(WMatrix) :: _CLASSNAME_(WMatrix)( const boost::python::object &_ob )
  {
    for( types::t_int i = 0; i < _DIM_; ++i )
    {
      boost::python::object a = _ob[i];
      for( types::t_int j = 0; j < _DIM_; ++j )
        x[i][j] = boost::python::extract<_TYPE_>( a[j] ); 
    }
  }

  std::ostream& operator<<( std::ostream &_o, const _CLASSNAME_(WMatrix)& _v )
  {
    _o << "( ";
    for( types::t_int i = 0; i < _DIM_; ++i )
    {
      _o << "[";
      for( types::t_int j = 0; j < _DIM_; ++j )
        _o << _v.x[i][j] << ( j < 2 ? ", ": "]" );
      _o << ( i <  _DIM_ - 1 ? ", ": " )" );
    }
    return _o;
  }
  std::string _CLASSNAME_(WMatrix) :: print() const
  {
    std::ostringstream sstr;
    sstr << *this;
    return sstr.str();
  }

  _CLASSNAME_(WMatrix)* operator*( const _CLASSNAME_(WMatrix) &_v,
                                  const _CLASSNAME_(WMatrix) &_b )
  {
    _CLASSNAME_(WMatrix) *result = new _CLASSNAME_(WMatrix);
    *( (FixedMatrix<_TYPE_, _DIM_>*) result) 
       = (FixedMatrix< _TYPE_, _DIM_ > ) _v * (FixedMatrix< _TYPE_, _DIM_ > ) _b;
    return result;
  }
  _CLASSNAME_(WMatrix)* operator*( const _CLASSNAME_(WMatrix) &_v, const _TYPE_ _b )
  {
    _CLASSNAME_(WMatrix) *result = new _CLASSNAME_(WMatrix);
    *( (FixedMatrix<_TYPE_, _DIM_>*) result) 
       = (FixedMatrix< _TYPE_, _DIM_ > ) _v * _b;
    return result;
  }
  _CLASSNAME_(WMatrix)* operator/( const _CLASSNAME_(WMatrix) &_v, const _TYPE_ _b )
  {
    _CLASSNAME_(WMatrix) *result = new _CLASSNAME_(WMatrix);
    *( (FixedMatrix<_TYPE_, _DIM_>*) result) 
       = (FixedMatrix< _TYPE_, _DIM_ > ) _v / _b;
    return result;
  }
  _CLASSNAME_(WMatrix)* operator*( const _TYPE_ _b, const _CLASSNAME_(WMatrix) &_v )
    { return operator*(_v, _b); }

  _CLASSNAME_(WVector)* operator*( const _CLASSNAME_(WMatrix) &_v,
                                  const _CLASSNAME_(WVector) &_b )
  {
    _CLASSNAME_(WVector) *result = new _CLASSNAME_(WVector);
    *( (FixedVector<_TYPE_, _DIM_>*) result) 
       = (FixedMatrix< _TYPE_, _DIM_ > ) _v * (FixedVector< _TYPE_, _DIM_ > ) _b;
    return result;
  }


}

#else

   class_< _CLASSNAME_(WVector) >( _PYTHONNAME_(Vector) )
     .def( init< object >() )
     .def( self += other<_CLASSNAME_(WVector)>() )
     .def( self + other<_CLASSNAME_(WVector)>() )
     .def( self - other<_CLASSNAME_(WVector)>() )
     .def( self -= other<_CLASSNAME_(WVector)>() )
     .def( other<_TYPE_>() * self )
     .def( self * other<_TYPE_>() )
     .def( self / other<_TYPE_>() )
     .def( self == other<_CLASSNAME_(WVector)>() )
     .def( other<_CLASSNAME_(WVector)>() == self )
     .def( self != other<_CLASSNAME_(WVector)>() )
     .def( other<_CLASSNAME_(WVector)>() != self )
     .def( "__str__", &_CLASSNAME_(WVector)::print );

   class_< _CLASSNAME_(WMatrix) >( _PYTHONNAME_(Matrix) )
     .def( init< object >() )
     .def( self += other<_CLASSNAME_(WMatrix)>() )
     .def( self -= other<_CLASSNAME_(WMatrix)>() )
     .def( self * other<_CLASSNAME_(WMatrix)>() )
     .def( self + other<_CLASSNAME_(WMatrix)>() )
     .def( self - other<_CLASSNAME_(WMatrix)>() )
     .def( self * other<_CLASSNAME_(WVector)>() )
     .def( other<_TYPE_>() * self )
     .def( self * other<_TYPE_>() )
     .def( self / other<_TYPE_>() )
     .def( self == other<_CLASSNAME_(WMatrix)>() )
     .def( other<_CLASSNAME_(WMatrix)>() == self )
     .def( self != other<_CLASSNAME_(WMatrix)>() )
     .def( other<_CLASSNAME_(WMatrix)>() != self )
     .def( "diag", &_CLASSNAME_(WMatrix)::diag )
     .def( "__str__", &_CLASSNAME_(WMatrix)::print );


#endif

#include "atat.py.hpp"

#endif 
