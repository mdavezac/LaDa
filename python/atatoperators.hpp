//
//  Version: $Id: atat.impl.hpp 897 2008-12-23 01:16:41Z davezac $
//

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <opt/types.h>

#define STREAM_VECTOR
#include <atat/fxvector.h>
#include <atat/vectmac.h>

namespace LaDa
{
  namespace atat
  {
    namespace details
    {
      template< class T_VECTOR >
        T_VECTOR add( const T_VECTOR& _a, const T_VECTOR& _b )
        {
          typedef typename Python::vector_introspection<T_VECTOR> :: type type;
          typedef atat::FixedVector<type, Python::vector_introspection<T_VECTOR> :: dim>  t_Fixed;
          return T_VECTOR( (const t_Fixed&)_a + (const t_Fixed&)_b ); 
        }
      template< class T_VECTOR >
        T_VECTOR minus( const T_VECTOR& _a, const T_VECTOR& _b )
        {
          typedef typename Python::vector_introspection<T_VECTOR> :: type type;
          typedef atat::FixedVector<type, Python::vector_introspection<T_VECTOR> :: dim>  t_Fixed;
          return T_VECTOR( (const t_Fixed&)_a - (const t_Fixed&)_b ); 
        }
      template< class T_VECTOR >
        T_VECTOR mul( const typename Python::vector_introspection<T_VECTOR> :: type _a,
                      const T_VECTOR& _b )
        {
          typedef typename Python::vector_introspection<T_VECTOR> :: type type;
          typedef atat::FixedVector<type, Python::vector_introspection<T_VECTOR> :: dim>  t_Fixed;
          return T_VECTOR( _a * (const t_Fixed&)_b ); 
        }

      template< class T_MATRIX >
        T_MATRIX addm( const T_MATRIX& _a, const T_MATRIX& _b )
        {
          typedef typename Python::matrix_introspection<T_MATRIX> :: type type;
          typedef atat::FixedMatrix<type, Python::matrix_introspection<T_MATRIX> :: dim>  t_Fixed;
          return T_MATRIX( (const t_Fixed&)_a + (const t_Fixed&)_b ); 
        }
      template< class T_MATRIX >
        T_MATRIX minusm( const T_MATRIX& _a, const T_MATRIX& _b )
        {
          typedef typename Python::matrix_introspection<T_MATRIX> :: type type;
          typedef atat::FixedMatrix<type, Python::matrix_introspection<T_MATRIX> :: dim>  t_Fixed;
          return T_MATRIX( (const t_Fixed&)_a - (const t_Fixed&)_b ); 
        }
      template< class T_MATRIX >
        T_MATRIX mulm( const typename Python::matrix_introspection<T_MATRIX> :: type _a,
                       const T_MATRIX& _b )
        {
          typedef typename Python::matrix_introspection<T_MATRIX> :: type type;
          typedef atat::FixedMatrix<type, Python::matrix_introspection<T_MATRIX> :: dim>  t_Fixed;
          return T_MATRIX( _a * (const t_Fixed&)_b ); 
        }
      template< class T_MATRIX >
        T_MATRIX mulm( const T_MATRIX& _a, const T_MATRIX& _b )
        {
          typedef typename Python::matrix_introspection<T_MATRIX> :: type type;
          typedef atat::FixedMatrix<type, Python::matrix_introspection<T_MATRIX> :: dim>  t_Fixed;
          return T_MATRIX( (const t_Fixed&)_a * (const t_Fixed&)_b ); 
        }
      template< class T_MATRIX >
        bool neq( const T_MATRIX& _a, const T_MATRIX& _b )
        {
          typedef typename Python::matrix_introspection<T_MATRIX> :: type type;
          typedef atat::FixedMatrix<type, Python::matrix_introspection<T_MATRIX> :: dim>  t_Fixed;
          return (const t_Fixed&)_a != (const t_Fixed&)_b; 
        }
      template< class T_MATRIX >
        typename Python::matrix_introspection< T_MATRIX > :: t_Vector 
          mulm( const T_MATRIX &_a, const typename Python::matrix_introspection<T_MATRIX> :: t_Vector &_b )
          {
            typedef typename Python::matrix_introspection<T_MATRIX> :: type type;
            typedef atat::FixedMatrix<type, Python::matrix_introspection<T_MATRIX> :: dim> t_Fixed;
            typedef atat::FixedVector<type, Python::matrix_introspection<T_MATRIX> :: dim> t_FixedVector;
            return (const t_Fixed&)_a * (const t_FixedVector&)_b; 
          }
      template< class T_MATRIX >
        typename Python::matrix_introspection< T_MATRIX > :: t_Vector 
          mulm( const typename Python::matrix_introspection<T_MATRIX> :: t_Vector &_b,
                const T_MATRIX &_a ) 
          {
            typedef typename Python::matrix_introspection<T_MATRIX> :: type type;
            typedef atat::FixedMatrix<type, Python::matrix_introspection<T_MATRIX> :: dim> t_Fixed;
            typedef atat::FixedVector<type, Python::matrix_introspection<T_MATRIX> :: dim> t_FixedVector;
            typename Python::matrix_introspection< T_MATRIX > :: t_Vector result(0,0,0);
            for(size_t i(0); i < 3; ++i) 
              for(size_t j(0); j < 3; ++j) 
                result(i) += _a(j,i) * _b(j);
            return result;
          }

      template< class T_MATRIX > T_MATRIX inv_rMatrix3d( const T_MATRIX &_a )
        { return !_a; }
      template< class T_MATRIX > T_MATRIX trans_rMatrix3d( const T_MATRIX &_a )
        { return ~_a; }
       
    } // namespace details.
                                      
#   define DECOPS( thistype ) \
      inline thistype operator+( const thistype& _a, const thistype& _b) \
        { return thistype( details::add( _a, _b ) ); } \
      inline thistype operator-( const thistype& _a, const thistype& _b) \
        { return thistype( details::minus( _a, _b ) ); } \
      inline thistype operator/( const thistype& _a,  \
                                 const Python::vector_introspection<thistype> :: type & _b) \
      { \
        typedef Python::vector_introspection<thistype> :: type type; \
        return thistype( details::mul( type(1) / _b, _a ) ); \
      } \
      inline thistype operator*( const thistype& _a, \
                                 const Python::vector_introspection<thistype> :: type & _b) \
        { return thistype( details::mul( _b, _a ) ); } \
      inline thistype operator*( const Python::vector_introspection<thistype> :: type & _a, \
                                 const thistype& _b ) \
        { return thistype( details::mul( _a, _b ) ); }

    DECOPS( rVector3d )
    DECOPS( iVector3d )

#   undef DECOPS

#   define DECOPS( thistype ) \
      inline thistype operator+( const thistype& _a, const thistype& _b) \
        { return thistype( details::addm( _a, _b ) ); } \
      inline thistype operator-( const thistype& _a, const thistype& _b) \
        { return thistype( details::minusm( _a, _b ) ); } \
      inline thistype operator/( const thistype& _a,  \
                                 const Python::matrix_introspection<thistype> :: type & _b) \
      { \
        typedef Python::matrix_introspection<thistype> :: type type; \
        return thistype( details::mulm( type(1) / _b, _a ) ); \
      } \
      inline thistype operator*( const thistype& _a, \
                                 const Python::matrix_introspection<thistype> :: type & _b) \
        { return thistype( details::mulm( _b, _a ) ); } \
      inline thistype operator*( const Python::matrix_introspection<thistype> :: type & _a, \
                                 const thistype& _b ) \
        { return thistype( details::mulm( _a, _b ) ); } \
      inline thistype operator*( const thistype & _a, const thistype& _b ) \
        { return thistype( details::mulm( _a, _b ) ); } \
      inline bool operator!=( const thistype & _a, const thistype& _b ) \
        { return details::neq( _a, _b ); }\
      inline Python::matrix_introspection<thistype> :: t_Vector \
        operator*( const thistype& _a,  const Python::matrix_introspection<thistype> :: t_Vector &_b) \
          { return details::mulm( _a, _b); } \
      inline Python::matrix_introspection<thistype> :: t_Vector \
        operator*( const Python::matrix_introspection<thistype> :: t_Vector &_b, \
                   const thistype& _a )  \
          { return details::mulm( _b, _a); }


    DECOPS( rMatrix3d )

#   undef DECOPS
  }
}

