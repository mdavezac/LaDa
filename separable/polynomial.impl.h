//
//  Version: $Id$
//

#if defined(__POLYFUNCTIONS__)

namespace  Separable
{
  namespace  Polynomial
  {
    namespace  detais
    {
#     if !defined( _NDEGREE_ )
#     define _NDEGREE_ 1 
      template< types::t_int D > struct Pow {};
      template<> struct Pow<0>
      {
        template< class T_TYPE >
          static T_TYPE pow( T_TYPE &_x ) { return T_TYPE(1); }
        template< class T_TYPE >
          static T_TYPE gradient( T_TYPE &_x ) { return T_TYPE(0); }
      };
      template<> struct Pow<1>
      {
        template< class T_TYPE >
          static T_TYPE pow( T_TYPE &_x ) { return _x; }
        template< class T_TYPE >
          static T_TYPE gradient( T_TYPE &_x ) { return T_TYPE(1); }
      };
      template<> struct Pow<-1>
      {
        template< class T_TYPE >
          static T_TYPE pow( T_TYPE &_x ) { return T_TYPE(1)/_x; }
        template< class T_TYPE >
          static T_TYPE gradient( T_TYPE &_x ) { return T_TYPE(-1)/Pow<2>::pow(_x); }
      };
#     elif _NDEGREE_ == 1
#     undef _NDEGREE_ 
#     define _NDEGREE_ 2
#     elif _NDEGREE_ == 2
#     undef _NDEGREE_ 
#     define _NDEGREE_ 3
#     elif _NDEGREE_ == 3
#     undef _NDEGREE_ 
#     define _NDEGREE_ 4
#     elif _NDEGREE_ == 4
#     undef _NDEGREE_ 
#     define _NDEGREE_ 5
#     elif _NDEGREE_ == 5
#     undef _NDEGREE_ 
#     define _NDEGREE_ 6
#     elif _NDEGREE_ == 6
#     undef _NDEGREE_ 
#     define _NDEGREE_ 7
#     elif _NDEGREE_ == 7
#     undef _NDEGREE_ 
#     define _NDEGREE_ 8
#     elif _NDEGREE_ == 8
#     undef _NDEGREE_ 
#     define _NDEGREE_ 9
#     elif _NDEGREE_ == 9
#     undef _NDEGREE_ 
#     define _NDEGREE_ 10
#     undef __POLYFUNCTIONS__
#     endif

  
#     if _NDEGREE_ > 1 
        template<> struct Pow<2*_NDEGREE_>
        {
          template< class T_TYPE >
            static T_TYPE pow( T_TYPE &_x )
              { T_TYPE a = Pow<_NDEGREE_>::pow(_x); return a*a; }
          template< class T_TYPE >
            static T_TYPE gradient( T_TYPE &_x )
              { return T_TYPE(2*_NDEGREE_)*Pow<2*_NDEGREE_-1>(_x);  }
        };
        template<> struct Pow<2*_NDEGREE_+1>
        {
          template< class T_TYPE >
            static T_TYPE pow( T_TYPE &_x ) { return Pow<2*_NDEGREE_>::pow(_x)*a; }
          template< class T_TYPE >
            static T_TYPE gradient( T_TYPE &_x )
              { return T_TYPE(2*_NDEGREE_+1)*Pow<2*_NDEGREE_>(_x);  }
        };
        template<> struct Pow<-2*_NDEGREE_>
        {
          template< class T_TYPE >
            static T_TYPE pow( T_TYPE &_x )
              { T_TYPE a = Pow<-1>::pow(_x); return Pow<2*_NDEGREE_>::pow(a); }
          template< class T_TYPE >
            static T_TYPE gradient( T_TYPE &_x )
              { return T_TYPE(-2*_NDEGREE_)/ (_x*Pow<2*_NDEGREE_>(_x));  }
        };
        template<> struct Pow<-2*_NDEGREE_+1>
        {
          template< class T_TYPE >
            static T_TYPE pow( T_TYPE &_x )
              { T_TYPE a = Pow<-1>::pow(_x); return Pow<2*_NDEGREE_+1>::pow(a); }
          template< class T_TYPE >
            static T_TYPE gradient( T_TYPE &_x )
              { return T_TYPE(-2*_NDEGREE_+1)/ (_x*Pow<2*_NDEGREE_+1>(_x));  }
        };
#     endif
    } // end of details namespace
  } // end of details Polynomial namespace
} // end of details Separable namespace

#endif // end of __POLYFUNCTIONS__
