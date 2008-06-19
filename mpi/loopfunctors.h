//
//  Version: $Id$
//
#ifndef _LOOP_FUNCTORS_H_
#define _LOOP_FUNCTORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace loop
{
  template< class T_OBJECT, class T_COMPONENT, class T_FUNCTOR>
  class  Component
  {
    public:
      typedef T_OBJECT t_Object;
      typedef T_COMPONENT t_Component;
      typedef T_FUNCTOR t_Functor;

    protected:
      t_Component *component;
      t_Functor functor;

    public:
      Component   ( t_Component *_comp, t_Functor _func )
                : component(_comp), functor(_func) {} 
      Component   ( const Component<t_Object, t_Component, t_Functor> &_comp )
                : component(_comp.component), functor(_comp.functor) {} 

      bool operator()( t_Object &_obj )
      { return _func( _obj.*component ); }
  };          

  template< class T_OBJECT,class T_FUNCTOR1, class T_FUNCTOR2 >
  class And
  {
    public:
      typedef T_OBJECT t_Object;
      typedef T_FUNCTOR1 t_Functor1;
      typedef T_FUNCTOR2 t_Functor2;

    protected:
      t_Functor1 functor1;
      t_Functor2 functor2;

    public:
      explicit
        And   ( t_Functor1 _func1, t_Functor2 _func2 )
            : functor1(_func1), functor2(_func2) {} 
      And   ( const And<t_Object, t_Functor1, t_Functor2> &_c )
          : functor1(_c.functor1), functor2(_c.functor2) {} 

      bool operator()( t_Object &_obj )
        { return _func1( _obj ) and _func2( _obj ); }
  };          

  template< class T_ITERATOR, class T_FUNCTOR>
  class OverIterators 
  {
    public:
      typedef T_ITERATOR t_Iterator;
      typedef T_FUNCTOR t_Functor;

    protected:
      t_Functor functor;

    public:
      OverIterators   ( t_Functor _func )
                    : functor(_func) {}
      OverIterators   ( const OverIterators< t_Iterator, t_Functor > &_c ) 
                    : functor(_c.functor) {}

      bool operator()( t_Iterator _first, t_Iterator _last )
      {
        bool result = true; 
        for (; _first != _last; ++_first ) 
          if ( not functor( *_first ) ) result = false;
        return false;
      }
  };          
  template< class T_ARRAY, t_int SIZE, class T_FUNCTOR>
  class OverFixedArray 
  {
    public:
      typedef T_ARRAY t_Array;
      typedef T_FUNCTOR t_Functor;

    protected:
      t_Functor functor;

    public:
      OverFixedArray   ( t_Functor _func )
                     : functor(_func) {}
      OverFixedArray   ( const OverFixedArray< t_Array, SIZE, t_Functor > &_c ) 
                     : functor(_c.functor) {}

      bool operator()( t_Array &_array )
      {
        bool result = true; 
        for (t_int i=0; i < SIZE; ++i ) 
          if ( not functor( _array[i] ) ) result = false;
        return false;
      }
  };          

  template< class T_CONTAINER, class T_FUNCTOR>
  class OverContainer 
  {
    public:
      typedef T_CONTAINER t_Container;
      typedef T_FUNCTOR t_Functor;

    protected:
      t_Functor functor;

    public:
      OverContainer   ( t_Functor _func )
                    : functor(_func) {}
      OverContainer   ( const OverContainer< t_Container, t_Functor > &_c ) 
                    : functor(_c.functor) {}

      bool operator()( t_Container &_obj )
      {
        typename t_Container :: iterator first = _obj.begin();
        typename t_Container :: iterator last = _obj.end();
        bool result = true;
        for (; first != last; ++first ) 
          if ( not functor( *first ) ) result = false;
        return false;
      }
  };          
  
} // namespace loop
//   typedef  mpi::broadcast_functor< types::t_real > t_broadcast
//   typedef OverFixedArray< types::t_real, 3, t_broadcast > t_OverFixedArray;
//   typedef Component< Crystal::Atom, atat::rVector3d, t_OverFixedArray > t_Component1;
//   typedef Component< Crystal::Structure::t_Atoms, t_broadcast >  t_Component2;
//   OverContainer< Crystal::Structure::t_Atoms, And< t_Component1, t_Component2 > >
//     serialize_structure( And( t_Component1( &Crystal::Atom::pos, 
//                                             OverFixedArray( t_broadcast( mpiobject ) ) ),
//                               t_Component2( &Crystal::Atom::type, t_broadcast( mpiobject ) ) ) );
//
//   serialize_structure( structure.atoms );
#endif
