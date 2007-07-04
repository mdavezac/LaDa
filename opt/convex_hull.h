#ifndef _CONVEX_HULL_BASE_H_
#define _CONVEX_HULL_BASE_H_

#include <list>
#include <functional>
#include <iostream>
#include <iomanip>

#include <tinyxml/tinyxml.h>

#include "opt/types.h"
#include "opt/opt_functors.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace ch
{
  template <class T_OBJECT>
  struct Vertex  // breaking point
  {
    typedef T_OBJECT t_Object;

    types::t_real y,x;
    t_Object object;

    Vertex() : y(0), x(0), object() {};
    template<class T_LOADOP>
    Vertex (const TiXmlElement &_node, T_LOADOP &_op) { Load( _node, _op ); }
    Vertex   (const types::t_real _y, const t_Object &_object) 
           : y(_y), x( _object.get_concentration() ), object(_object) {}
    Vertex(const Vertex<t_Object> &_v ) : y(_v.y), x(_v.x), object(_v.object) {};

    bool vx_less(const Vertex<t_Object> &_v) const
      { return (x < _v.x); }
    bool x_less(const Vertex<t_Object> &_v) const
      { return (x < _v.x); }
    bool x_less(const types::t_real _x) const
      { return (x < _x); }
    bool x_greater(const types::t_real _x) const
      { return (x > _x); }
    bool x_geq(const types::t_real _x) const // greater or equal
      { return ( std::abs(x-_x) < types::tolerance ) ? true : (x > _x); }
    bool x_less(const Vertex<t_Object> &_v1, const Vertex<t_Object> &_v2) const
      { return (_v1.x < _v2.x); }
    bool E_less(const Vertex<t_Object> &_v) const
      { return (y < _v.y); }
    bool E_less(const Vertex<t_Object> &_v1, const Vertex<t_Object> &_v2) const
      { return (_v1.y < _v2.y); }
    bool operator==(const Vertex<t_Object> &_v) const
    {
      if (    std::abs(y - _v.y) > types::tolerance
           or std::abs(x - _v.x ) > types::tolerance )
        return false;
      return (object == _v.object); 
    }

    void print_out(std::ostream &stream) const
      { stream << x << " " << y << std::endl; }
    template< class T_LOADOP >
    bool Load( const TiXmlElement &_node, T_LOADOP &_op );
    template< class T_SAVEOP >
    void Save( TiXmlElement &_node, T_SAVEOP &_op ) const;

    types::t_real get_concentration() const
      { return x; }
  };

  struct HalfLine 
  { 
    types::t_real a, b;
    types::t_real x_end;
    types::t_real evaluate(const types::t_real x) const { return (a*x + b); }
    types::t_real evaluate_gradient() const { return a; }
    HalfLine() : a(0), b(0), x_end(0) {};
    HalfLine( const HalfLine &_s ) : a(_s.a), b(_s.b), x_end(_s.x_end) {};
    template<class T_OBJECT>
      HalfLine( const Vertex<T_OBJECT> &_s, const Vertex<T_OBJECT> &_e)
      {
        b = _s.x - _e.x; 
        if ( b == 0 ) 
          { std::cerr << "Error initialising linear segment" << std::endl; exit(0); }
        a = (_s.y - _e.y)/b; b = _e.y - _e.x*a;
        x_end = _e.x;
      }
    bool x_greater(const types::t_real _x) const
      { return (_x < x_end); }
  };

  template<class T_OBJECT>
  class Base
  {
    public:
      typedef T_OBJECT t_Object;
      typedef Vertex<t_Object> t_Vertex;
      typedef std::list< Vertex<t_Object> > t_Vertices;
      typedef std::list< HalfLine > t_HalfLines;

    protected:
      t_Vertices vertices;
      t_HalfLines segments;

    public:
      Base() {}
      Base   ( const Base<t_Object>& _b) 
           : vertices(_b.vertices), segments(_b.segments) {}
      ~Base() {}

      bool add( const types::t_real _y, const t_Object &_o, bool _do_check = true );

      types::t_real evaluate(const types::t_real _x) const;
      types::t_real evaluate_gradient( const types::t_real _x ) const;
      void clear() { vertices.clear(); }
      typename std::list< t_Vertex > :: const_iterator begin_vertex() const
        { return vertices.begin(); }
      typename std::list< t_Vertex > :: const_iterator end_vertex() const
        { return vertices.end(); }
      template< class T_LOADOP >
      bool Load(const TiXmlElement &_node, T_LOADOP &_op);
      template< class T_SAVEOP >
      void Save( TiXmlElement &_node, T_SAVEOP &_op ) const;

#ifdef _MPI
      bool broadcast( mpi::BroadCast &_bc )
      {
        types::t_int n = vertices.size();
        if ( not _bc.serialize(n) ) return false;
        if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
          vertices.resize(n);
        typename t_Vertices :: iterator i_vertex = vertices.begin();
        typename t_Vertices :: iterator i_vertex_end = vertices.end();
        for(; i_vertex != i_vertex_end; ++i_vertex )
          if(    ( not _bc.serialize<t_Object>( i_vertex->object ) ) 
              or ( not _bc.serialize( i_vertex->x ) )
              or ( not _bc.serialize( i_vertex->y ) )       ) return false;

        if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
        {
          while( not weed_out() ) {}
          build_function();
        }

        return true;
      }
#endif

    protected:
      void build_function();
      bool weed_out(); // removes unconvex points one by one
  };
      

  template<class T_OBJECT>
  bool Base<T_OBJECT> :: add( const types::t_real _y, const t_Object &_o, bool _do_check )
  {
     t_Vertex vertex(_y,_o);
     // empty vertex list!!
     if ( vertices.empty() )
     {
       vertices.push_back(vertex);
       return true;
     }

     typename t_Vertices :: iterator i_begin = vertices.begin();
     typename t_Vertices :: iterator i_end = vertices.end();
     typename t_Vertices :: iterator i_insert = i_begin; ++i_insert;

     // exception case where there is only one point on 
     // the convex hull yet
     if ( i_insert == i_end )
     {
       if ( std::abs(i_begin->x - vertex.x) < types::tolerance )
       { 
         if ( i_begin->y > vertex.y )
           *i_begin = vertex;
         return false; // only one point, no convexhull yet!
       }
       else 
         i_begin->x > vertex.x ? vertices.push_front( vertex ) : vertices.push_back( vertex );
       build_function();
       return true;
     }

     std::const_mem_fun1_ref_t<bool, t_Vertex, const types::t_real> func(&t_Vertex::x_geq );
     
     // does not change convex hull if this is true
     if ( vertex.y - evaluate(vertex.x) > 0.0 )
       return false; 
     
     // first finds where to insert new vertex
     i_insert = std::find_if( i_begin, i_end, std::bind2nd(func, vertex.x) );
     if ( i_insert == i_end ) // outside of range
       vertices.push_back(vertex);
     else if( std::abs(i_insert->x - vertex.x) < types::tolerance ) 
     {  // point already exists at same concentration
       if ( std::abs( i_insert->y - vertex.y ) < types::tolerance ) return false;
       *i_insert = vertex;
     }
     else // otherwise, inserts new vertex
       vertices.insert(i_insert, vertex );
     
     // makes convex-hull convex
     while ( not ( _do_check and weed_out() ) );
     
     // rebuilds our step-wise function
     build_function();
     
     return true;
  }
  
  // builds linear segment convex_hull
  template<class T_OBJECT>
  void Base<T_OBJECT> :: build_function ()
  {
    if ( vertices.size() < 2 )
      return;

    segments.clear(); // builds from scratch
      
    typename t_Vertices :: iterator i_v = vertices.begin();
    typename t_Vertices :: iterator i_v2 = i_v;
    typename t_Vertices :: iterator i_end = vertices.end();

    for ( ++i_v; i_v != i_end; ++i_v, ++i_v2 )
      segments.push_back( HalfLine(*(i_v2), *i_v) );
  }

  template<class T_OBJECT>
  bool Base<T_OBJECT> :: weed_out()
  {
    typename t_Vertices :: iterator i_begin, i_end, i_mid, i_last;
    types::t_real p1, p2;

    if (vertices.size() < 3 )  // two point hull always convex
      return true;

    // the following loop checks that successive triangles are convex
    i_begin = vertices.begin();
    i_last = vertices.end();
    i_mid = i_begin; ++i_mid;
    i_end = i_mid; ++i_end;
    while ( i_end != i_last ) 
    {
      p1 = ( i_begin->y - i_mid->y ) * ( i_mid->x - i_end->x );
      p2 = ( i_mid->y - i_end->y ) * ( i_begin->x - i_mid->x );
      if ( p1 >= p2 )
        break;
      ++i_begin; ++i_mid; ++i_end;
    }
    if ( i_end == i_last )  // all triangles are convex
      return true;

    // one triangle is concave -> erase midpoint
    vertices.erase(i_mid);
    
    // we don't know wether this is convex yet -> returns false
    return false; 
  }

  template<class T_OBJECT> template <class T_LOADOP>
  bool Base<T_OBJECT> :: Load( const TiXmlElement &_node, T_LOADOP &_op )
  {
    std::string name = _node.Value();
    const TiXmlElement *child;
    if ( name.compare("ConvexHull" ) != 0 )
      child = _node.FirstChildElement("ConvexHull");
    else
      child = &_node;
    if ( not child )
      return false;

    child = child->FirstChildElement( "Vertex" );
    std::vector< t_Vertex > sortcopy;
    for ( ; child; child=child->NextSiblingElement( "Vertex" ) )
    {
      t_Vertex vertex(*child, _op);
      sortcopy.push_back( vertex );
    }
    switch ( sortcopy.size() )
    {
      case 0: return false; break;
      case 1: vertices.push_back( sortcopy.back() ); return true; break;
      default: break;        
    }

    std::sort( sortcopy.begin(), sortcopy.end(), std::mem_fun_ref( &Vertex<t_Object>::vx_less ) );
    vertices.resize( sortcopy.size() );
    std::copy( sortcopy.begin(), sortcopy.end(), vertices.begin() );


    while( !weed_out() );

    build_function();

    return true;
  }

  template<class T_OBJECT> template< class T_LOADOP >
  bool Vertex<T_OBJECT> :: Load( const TiXmlElement &_node, T_LOADOP &_op )
  {
    std::string name = _node.Value();
    if ( name.compare("Vertex" ) != 0 )
      return false;

    _node.Attribute("y", &y);
    _node.Attribute("x", &x);
    return _op(object, _node);
  }

  template<class T_OBJECT> template<class T_SAVEOP>
  void Base<T_OBJECT> :: Save( TiXmlElement &_node, T_SAVEOP &_op ) const
  {
    TiXmlElement *parent;
    
    parent = new TiXmlElement("ConvexHull");
    _node.LinkEndChild(parent);

    typename t_Vertices :: const_iterator i_v = vertices.begin();
    typename t_Vertices :: const_iterator i_end = vertices.end();
    for( ; i_v != i_end; ++i_v )
      i_v->Save(*parent, _op);
  }
  template<class T_OBJECT> template<class T_SAVEOP>
  void Vertex<T_OBJECT> :: Save( TiXmlElement &_node, T_SAVEOP &_op ) const
  {
    TiXmlElement *child = new TiXmlElement("Vertex");
    child->SetDoubleAttribute( "y", y );
    child->SetDoubleAttribute( "x", x );
    _op( object, *child );
    _node.LinkEndChild(child);
  }

  template<class T_OBJECT>
  types::t_real Base<T_OBJECT> :: evaluate( const types::t_real _x) const
  {
    if ( segments.size() == 0 )
      return 0.0;

    t_HalfLines :: const_iterator i_begin = segments.begin();
    t_HalfLines :: const_iterator i_end = segments.end();
    t_HalfLines :: const_iterator i_which;
    std::const_mem_fun1_ref_t<bool, HalfLine, const types::t_real> 
                             func( &HalfLine :: x_greater );
    i_which = std::find_if( i_begin, i_end, std::bind2nd( func, _x ) );

    if ( i_which == i_end ) // exception case ...
      --i_which;

    return i_which->evaluate(_x);
  }

  template<class T_OBJECT>
  types::t_real Base<T_OBJECT> :: evaluate_gradient( const types::t_real _x) const
  {
    if ( segments.empty() )
      return 0.0;

    typename t_HalfLines :: const_iterator i_begin = segments.begin();
    typename t_HalfLines :: const_iterator i_end = segments.end();
    typename t_HalfLines :: const_iterator i_which;
    std::const_mem_fun1_ref_t<bool, HalfLine, const types::t_real> 
                             func( &HalfLine :: x_greater );
    i_which = std::find_if( i_begin, i_end, std::bind2nd( func, _x ) );

    if ( i_which == i_end ) // exception case ...
      --i_which; 

    return i_which->a;
  }

} // namespace ch

#endif //  _CONVEX_HULL_BASE_H_
