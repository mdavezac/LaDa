//
//  Version: $Id$
//

namespace LaDa
{
  namespace opt
  {

    namespace ConvexHull
    {
      template<class T_OBJECT> template<class ARCHIVE>
        void Base<T_OBJECT> :: save(ARCHIVE & _ar, const unsigned int _version) const
        {
          _ar & vertices;
          while( not weed_out() ) {}
          build_function();
        
          return true;
        }

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
         typename t_Vertices :: iterator i_insert;

         // exception case where there is only one point on 
         // the convex hull yet
         if ( vertices.size() == 1 )
         {
           if ( Fuzzy::eq(i_begin->x, vertex.x) )
           { 
             if ( Fuzzy::gt(i_begin->x, vertex.x) )
               *i_begin = vertex;
             return false; // only one point, no convexhull yet!
           }
           else if ( i_begin->x > vertex.x ) vertices.push_front( vertex );
           else vertices.push_back( vertex );
           build_function();
           return true;
         }

         std::const_mem_fun1_ref_t<bool, t_Vertex,
                                   const types::t_real> func(&t_Vertex::x_geq );
         
         // does not change convex hull if this is true
         if ( Fuzzy::gt(vertex.y, evaluate(vertex.x) ) )
           return false; 
         
         // first finds where to insert new vertex
         i_insert = std::find_if( i_begin, i_end, std::bind2nd(func, vertex.x) );
         if ( i_insert == i_end ) // outside of range
           vertices.push_back(vertex);
         else if( Fuzzy::eq(i_insert->x, vertex.x ) )
         {  // point already exists at same concentration
           if (  Fuzzy::eq(i_insert->y, vertex.y ) ) return false;
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
      bool Base<T_OBJECT> :: Load( const TiXmlElement &_node, const T_LOADOP &_op )
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

        std::sort( sortcopy.begin(), sortcopy.end(), std::mem_fun_ref( &Vertex<t_Object>::x_less ) );
        vertices.resize( sortcopy.size() );
        std::copy( sortcopy.begin(), sortcopy.end(), vertices.begin() );


        while( !weed_out() );

        build_function();

        return true;
      }

      template<class T_OBJECT> template< class T_LOADOP >
      bool Vertex<T_OBJECT> :: Load( const TiXmlElement &_node, const T_LOADOP &_op )
      {
        std::string name = _node.Value();
        if ( name.compare("Vertex" ) != 0 )
          return false;

        _node.Attribute("y", &y);
        _node.Attribute("x", &x);
        return _op(object, _node);
      }

      template<class T_OBJECT> template<class T_SAVEOP>
      bool Base<T_OBJECT> :: Save( TiXmlElement &_node, const T_SAVEOP &_op ) const
      {
        TiXmlElement *parent;
        
        parent = new TiXmlElement("ConvexHull");
        if( not parent  ) return false;
        _node.LinkEndChild(parent);

        typename t_Vertices :: const_iterator i_v = vertices.begin();
        typename t_Vertices :: const_iterator i_end = vertices.end();
        for( ; i_v != i_end; ++i_v )
          if ( not i_v->Save(*parent, _op) ) break;

        return i_v == i_end;
      }
      template<class T_OBJECT> template<class T_SAVEOP>
      bool Vertex<T_OBJECT> :: Save( TiXmlElement &_node, const T_SAVEOP &_op ) const
      {
        TiXmlElement *child = new TiXmlElement("Vertex");
        if ( not child ) return false;
        child->SetDoubleAttribute( "y", y );
        child->SetDoubleAttribute( "x", x );
        _op( object, *child );
        _node.LinkEndChild(child);
        return true;
      }

      template<class T_OBJECT>
      types::t_real Base<T_OBJECT> :: evaluate( const types::t_real _x) const
      {
        if ( segments.size() == 0 ) return 0.0;

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
        if ( segments.empty() ) return 0.0;

        typename t_HalfLines :: const_iterator i_end = segments.end();
        typename t_HalfLines :: const_iterator i_which;
        std::const_mem_fun1_ref_t<bool, HalfLine, const types::t_real> 
                                 func( &HalfLine :: x_greatereq );
        i_which = std::find_if( segments.begin(), i_end, std::bind2nd( func, _x ) );

        // exception case ...
        if ( i_which == i_end ) --i_which; 
        else if ( Fuzzy::eq( i_which->x_end, _x ) )
          return 0e0;

        return i_which->a;
      }

      template<class T_OBJECT>
      types::t_real Base<T_OBJECT> :: evaluate_left_gradient( const types::t_real _x) const
      {
        if ( segments.empty() ) return 0.0;

        typename t_HalfLines :: const_iterator i_end = segments.end();
        typename t_HalfLines :: const_iterator i_which;
        std::const_mem_fun1_ref_t<bool, HalfLine, const types::t_real> 
                                 func( &HalfLine :: x_greatereq );
        i_which = std::find_if( segments.begin(), i_end, std::bind2nd( func, _x ) );

        // exception case ...
        if ( i_which == i_end ) --i_which;

        return i_which->a;
      }

      template<class T_OBJECT>
      types::t_real Base<T_OBJECT> :: evaluate_right_gradient( const types::t_real _x) const
      {
        if ( segments.empty() ) return 0.0;

        typename t_HalfLines :: const_iterator i_end = segments.end();
        typename t_HalfLines :: const_iterator i_which;
        std::const_mem_fun1_ref_t<bool, HalfLine, const types::t_real> 
                                 func( &HalfLine :: x_greatereq );
        i_which = std::find_if( segments.begin(), i_end, std::bind2nd( func, _x ) );

        // exception case ...
        if ( i_which == i_end ) --i_which;
        else
        {
          if ( Fuzzy::eq( i_which->x_end, _x ) ) ++i_which;
          if ( i_which == i_end ) --i_which;
        }

        return i_which->a;
      }

    } // namespace ConvexHull
  } // namespace opt
} // namespace LaDa
