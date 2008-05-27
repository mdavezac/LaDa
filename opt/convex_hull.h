//
//  Version: $Id$
//
#ifndef _CONVEX_HULL_BASE_H_
#define _CONVEX_HULL_BASE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <list>
#include <functional>
#include <iostream>
#include <iomanip>

#include <tinyxml/tinyxml.h>

#include <mpi/mpi_object.h>

#include "types.h"
#include "function_functors.h"
#include "fuzzy.h"

#include <mpi/mpi_object.h>

#include <boost/serialization/serialization.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <print/stdout.h>

//! trash-can namespace for anything that doesn't quite go anywhere
namespace opt
{

/**  \brief Contains template class which can will make a convex-hull out of
            any object of type \a T_OBJECT.
     \details First, let us define  what is a convex-hull.
      Imagine a space comprised of \a T_OBJECTs \f$\sigma\f$ to each of which
      can be attached an "energy" \f$q_\sigma\f$ and a
      concentration \f$x_\sigma\f$ in the range \f$x_\sigma \in [0,1]\f$. A
      convex-hull can be constructed in two steps.  First, one selects at each
      concentration \f$x_i\f$ the \a T_OBJECT \f$\sigma_i\f$ with lowest
      energy \f$q_\sigma\f$. Then, one deletes any \a T_OBJECT
      \f$\sigma_i\f$ that can disproportionate into a sum of two
      neighboring \a T_OBJECT \f$\sigma_{i-1}\f$ and
      \f$\sigma_{i+1}\f$, with
      \f[x_{\sigma_{i-1}} < x_{\sigma_i} < x_{\sigma_{i+1}}\f], 
      and 
      \f[
          q_{\sigma_i} > \frac{x_{\sigma_{i+1}} -
          x_{\sigma_i}}{x_{\sigma_{i+1}}-x_{\sigma_{i-x}}}
          q_{\sigma_{i-1}}\\+ \frac{x_\sigma -
          x_{\sigma_{i-1}}}{x_{\sigma_{i+1}}-x_{\sigma_{i-1}}} q_{\sigma_{i+1}} 
      \f].
      The \a T_OBJECT which satisfy the above conditions are known as
      the breaking-points of the convex-hull. Using these breaking-points, we
      can construct a piece-wise linear function \f$C(\sigma)\f$ connecting the
      breaking points in the \f$(q,x)\f$ plane. \f$C(\sigma)\f$ is the convex-hull.

      The class opt::ConvexHull::Base is an implementation of the convex-hull
      algorithm described above. More specifically, when fed an \a T_OBJECT A,
      it checks it is a breaking-point with respect to previously fed \a
      T_OBJECTS. If A is indeed a breaking-point, opt::ConvexHull::Base
      stores it for future reference. As can be seen, obaining the true
      convex-hull of a space requires, in theory, that you feed each every
      object in the space to opt::ConvexHull::Base. In practice, one could use
      some sampling method...

      ConvexHull::Base defines a function which takes a scalar on input. It is
      piecewise continuous. As a result its left and right derivatives exist
      for all points. For convenience a \e derivative is also defined which
      returns zero at the positions of breakpoints.

      The <em>x</em>-axis is determined by \a T_OBJECT::get_concentration(),
      The <em>y</em>-axis should be given on input when submitting an instance
      of \a T_OBJECT to a ConvexHull::Base object, <em> via</em>
      ConvexHull::Base::add(). Once constructed a convex-hull is nothing more
      than a single variable function. 
     
       \par Required behaviors for \a T_OBJECT
       The following behaviors are required from \a T_OBJECT and should be
      implemented to use it in ConvexHull,
       - get_concentration() returns a types::real which should represent the concentration
       - copy constructor or operator=(T_OBJECT&, const TOBJECT&) 
       - a functor capable of saving a \a T_OBJECT  to an XML node, eg GA::SaveObject 
       - a functor capable of loading a \a T_OBJECT  from an XML node, eg GA::LoadObject
       . 
      
       \par Design Philosophy
       Convexhulls are defined using two template classes and a single non-template class:
         - Vertex defines  a breakpoint in the convex hull through an copied
           instance of T_OBJECT, ConvexHull::Vertex::object, and two stored variables,
           ConvexHull::Vertex::x and ! ConvexHull::Vertex::y 
         - HalfLine defines a half-line in the \f$(x,y)\f$ plane. 
         - Base handles a collection of vertices. It checks whether the instances
           of \a T_OBJECT are breaking points and should be added to the ConvexHull::Vertex
           collection. The vertices are used to create a collection of
           ConvexHull::HalfLine objects which define the convex-hull in the \f$(x,y)\f$
           plane.
         . 
         Once created, a convex hull can be evaluated <em>via</em>
         ConvexHull::Base::evaluate(). It can be loaded and saved to XML
         <em>via</em> ConvexHull::Base::Load() and ConvexHull::Base::Save().
     
      Here is a sample code. The full implementation is not given.
      The header file of the T_OBJECT could go as follows,
      \code
        // The T_OBJECT type that will be stored as a convex-hull
        class Object
        {
            types::t_real some_member;
            std::vector<types::t_real> X;
          public: 
            //! Copy constructor
            Object(const Object& _c ) : some_member(_c.some_member), X(_c.X) {}
            //! returns object concentration  
            types::t_real get_concentration() const
            {
             return std::accumulate( X.begin(), X.end(), 0.0 ) / (types::t_real) X.size();
            }
        };
     
        //! Loads an Object from XML
        class LoadFunctor
        {
          public:
          bool operator()(Object&, const TiXmlElement &_node); 
        };
        //! Saves an Object to XML
        class SaveFunctor
        {
          public:
          bool operator()(const Object&, const TiXmlElement &_node); 
        };
      \endcode
      The convex hull itself, some object, and this object's y coordinate are declared
      \code
         ConvexHull::Base<Object> convexhull;
         Object object;
         types::t_real y;
      \endcode
      Something is done with \a object and with \a y, and then \a object can be
      added to the convex-hull,
      \code
         convexhull.add( y, object); // submits object to convexhull
         convexhull.evaluate(0.5);  // evaluates convexhull 
      \endcode

      \xmlrestart The convex-hull saves and reloads itself to and from XML in
      the following format:
      \code
         <ConvexHull>
            <Vertex x="?" y="?"/>
              ... T_OBJECT as saved by template parameter T_SAVEOP &_op of ConvexHull::Base::Save
            </Vertex
            <Vertex x="?" y="?"/>
              ... T_OBJECT as saved by template parameter T_SAVEOP &_op of ConvexHull::Base::Save
            </Vertex
            ... other vertices
         </ConvexHull>
      \endcode The attributes \a x and \a y of the \<%Vertex\> store the
      concentration (\a x) and the "energy" (\a y) of that Vertex. 
 */
namespace ConvexHull
{
  //! \brief Represents a Breaking-Point in the convex-hull
  //! \details It holds an instanciation of a \a T_OBJECT in Vertex::Object, as
  //! well as the coordinates of Vertex::x and Vertex::y of that breaking-point
  template <class T_OBJECT>
  struct Vertex  
  {
    typedef T_OBJECT t_Object; //!< the object type

    types::t_real x; //!< the \a x coordinate in the convex-hull plane
    types::t_real y; //!< the \a y coordinate in the convex-hull place
    t_Object object; //!< a copied instance of a t_Object

    Vertex() : y(0), x(0), object() {}; //!< Constructor
    //! \brief Constructor capable of loading itself from an XML node
    //! \param _node XML node from which to load the Vertex
    //! \param _op Functor capable of loading an instance of t_Object
    template<class T_LOADOP>
    Vertex (const TiXmlElement &_node, T_LOADOP &_op) { Load( _node, _op ); }
    //! \brief Constructor and Initializer
    //! \details Vertex::x is set using \a _object.get_concentration()
    Vertex   (const types::t_real _y, const t_Object &_object) 
           : y(_y), x( _object.get_concentration() ), object(_object) {}
    //! copy constructor
    Vertex(const Vertex<t_Object> &_v ) : y(_v.y), x(_v.x), object(_v.object) {};

    //! defines an operator< for the Vertex::x coordinate
    bool x_less(const Vertex<t_Object> &_v) const
      { return (x < _v.x); }
    //! defines an operator> for the Vertex::x coordinate
    bool x_geq(const types::t_real _x) const // greater or equal
      { return ( std::abs(x-_x) < types::tolerance ) ? true : (x > _x); }
    //! defines an operator< for the Vertex::y coordinate
    bool E_less(const Vertex<t_Object> &_v) const
      { return (y < _v.y); }
    //! defines an operate== which compares Vertex::x, Vertex::y, and Vertex::object
    bool operator==(const Vertex<t_Object> &_v) const
    {
      if (    std::abs(y - _v.y) > types::tolerance
           or std::abs(x - _v.x ) > types::tolerance )
        return false;
      return (object == _v.object); 
    }

    //! prints out Vertex::x and Vertex::y on single line to \a _stream.
    void print_out(std::ostream &_stream) const
      { _stream << x << " " << y << std::endl; }
    //! \brief Loads a Vertex from XML
    //! \param _node  XML node from which to load the Vertex
    //! \param _op a functor capable of loading a Vertex::t_Object from XML
    template< class T_LOADOP >
    bool Load( const TiXmlElement &_node, T_LOADOP &_op );
    //! \brief Save to XML
    //! \param _node  XML node to which to save the Vertex
    //! \param _op a functor capable of saving a Vertex::t_Object to XML
    template< class T_SAVEOP >
    bool Save( TiXmlElement &_node, T_SAVEOP &_op ) const;

    //! return Vertex::x
    types::t_real get_concentration() const
      { return x; }
    //! Serializes a vertex.
    template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
      { _ar & x; _ar & y; _ar & object; }
  };
  //! Prints out vertex using usual << operator.
  template< class T_OBJECT >
  std::ostream& operator<<( std::ostream &_os, const Vertex<T_OBJECT> &_v )
    { _v.print_out( _os ); return _os; }

  //! \brief Defines a halfline in the \f$(x,y)\f$ plane
  //! \details The half line is defined using three values,
  //!  - Vertex::a is the slope of the line
  //!  - Vertex::b is an offset along the <em>y</em>-axis
  //!  - Vertext::x_end defines the end-point of the half-line. The halfline exists for
  //! \f$x < x_\mathrm{end}\f$
  struct HalfLine 
  { 
    types::t_real a; //!< Slope
    types::t_real b; //!< <em>y</em>-axis offset
    types::t_real x_end; //!< rightmost end-point
    //! Returns the value of the line (even if beyond end-point) for \a x
    types::t_real evaluate(const types::t_real x) const { return (a*x + b); }
    //! Returns the gradient, eg the slope HalfLine::a
    types::t_real evaluate_gradient() const { return a; }
    //! Constructor
    HalfLine() : a(0), b(0), x_end(0) {};
    //! Copy Constructor
    HalfLine( const HalfLine &_s ) : a(_s.a), b(_s.b), x_end(_s.x_end) {};
    //! \brief Constructor using two Vertex objects
    //! \details the halfline connects \a _s and \a _e, with \a _e the endpoint
    template<class T_OBJECT>
      HalfLine( const Vertex<T_OBJECT> &_s, const Vertex<T_OBJECT> &_e)
      {
        b = _s.x - _e.x; 
        if ( b == 0 ) 
          { std::cerr << "Error initialising linear segment" << std::endl; exit(0); }
        a = (_s.y - _e.y)/b; b = _e.y - _e.x*a;
        x_end = _e.x;
      }
    //! returns \a _x < Vertex::x_end 
    bool x_greater(const types::t_real _x) const
      { return _x < x_end; }
    //! returns \a _x <= Vertex::x_end 
    bool x_greatereq(const types::t_real _x) const
      { return _x < x_end or Fuzzy::eq( _x, x_end); }
  };

  //! \brief Handles a collection of Vertex objects defining the convex-hull
  //! \details It handles addition and removal of submitted objects, such that
  //! the hull is kept convex and as low in \a y as possible.
  //! The evaluation of the convex-hull is done using a collection of HalfLine
  //! objects constructed from the vertices. 
  template<class T_OBJECT>
  class Base
  {
    public:
      typedef T_OBJECT t_Object; //!< Object type
      typedef Vertex<t_Object> t_Vertex;//!< Vertex type
      //! \brief Vertex collection type
      //! \details t_Vertices is an std::list since adding and removing
      //! vertices should happen relatively often as the convex-hull is refined.
      //! \warning Technically, one should not expect orderto preserved in
      //! lists. Nonetheless, this implementation does... The advantage of
      //! using a list is that we can insert objects in amortized time.
      typedef std::list< Vertex<t_Object> > t_Vertices; 
      //! \brief HalfLine collection type
      //! \details Just as t_Vertices, this is a list since half-lines are
      //! expected to be added and removed relatively often as the convex-hull is
      //! refined. 
      //! \warning Technically, one should not expect orderto preserved in
      //! lists. Nonetheless, this implementation does... The advantage of
      //! using a list is that we can insert objects in amortized time.
      typedef std::list< HalfLine > t_HalfLines;

    protected:
      t_Vertices vertices; //!< Vertex collection
      t_HalfLines segments; //!< HalfLine collection

    public:
      Base() {} //!< Constuctor
      //! Copy Constructor
      Base   ( const Base<t_Object>& _b)
           : vertices(_b.vertices), segments(_b.segments) {}
      ~Base() {} //!< Destructor

      //! \brief Tries and adds a t_Object \a _o to the vertices
      //! \details \a _o is added only if it is lower in energy than the convex-hull.
      //! If it is added, breaking-point which do not fulfill the convexity
      //! requirements are removed.
      //! \param _y energy of the t_Object \a _o
      //! \param _o t_Object to try and add to the convex-hull
      //! \param _do_check whether to check for convexity if and after \a _o has been added
      bool add( const types::t_real _y, const t_Object &_o, bool _do_check = true );

      //! return the value of the convex-hull at \a _x
      types::t_real evaluate(const types::t_real _x) const;
      //! \brief return the gradient of the convex-hull at \a _x
      //! \details Returns 0 if \e _x is at a breaking point
      types::t_real evaluate_gradient( const types::t_real _x ) const;
      //! \brief return the gradient of the convex-hull at \a _x
      //! \details At a breaking point, returns the left derivative.
      types::t_real evaluate_left_gradient( const types::t_real _x ) const;
      //! \brief return the gradient of the convex-hull at \a _x
      //! \details At a breaking point, returns the left derivative.
      types::t_real evaluate_right_gradient( const types::t_real _x ) const;
      //! removes all vertices
      void clear() { vertices.clear(); }
      //! returns an iterator to the beginning of the Vertex collection
      typename std::list< t_Vertex > :: const_iterator begin_vertex() const
        { return vertices.begin(); }
      //! returns an iterator to the end of the Vertex collection
      typename std::list< t_Vertex > :: const_iterator end_vertex() const
        { return vertices.end(); }
      //! \brief Loads the convex-hull from XML
      //! \param _node  XML node from which to load the Vertex
      //! \param _op a functor capable of loading a t_Object from XML
      template< class T_LOADOP >
      bool Load(const TiXmlElement &_node, T_LOADOP &_op);
      //! \brief Saves convex-hull to XML
      //! \param _node  XML node to which to save the Vertex
      //! \param _op a functor capable of saving a t_Object to XML
      template< class T_SAVEOP >
      bool Save( TiXmlElement &_node, T_SAVEOP &_op ) const;

      //! Returns number of breakpoints in convex-hull
      types::t_unsigned size() const { return vertices.size(); }

      //! loads a convex-hull.
      template<class ARCHIVE> void load(ARCHIVE & _ar, const unsigned int _version)
        { _ar & vertices; }
      //! saves a convex-hull.
      template<class ARCHIVE> void save(ARCHIVE & _ar, const unsigned int _version) const;
      BOOST_SERIALIZATION_SPLIT_MEMBER()
      //! Returns  a string containing the convex-hull in xmgrace .agr format.
      std::string print() const
      { 
        std::ostringstream sstr;
        std::for_each( vertices.begin(), vertices.end(),
                       boost::lambda::var(sstr) 
                         << boost::lambda::bind( &Vertex<T_OBJECT>::object, boost::lambda::_1 )
                         << boost::lambda::constant( " " )
                         << boost::lambda::_1  );
        sstr << "&";
        return sstr.str();
      }

    protected:
      //! Builds a collection of HalfLine from the vertices
      void build_function();
      //! \brief Removes the first break-point it finds which makes the convex-hull not convex
      //! \details If a break-point is removed, then it returns false. If no
      //! unconvex breaking-point were found, than it returns true. This
      //! function should be called over and over again until it returns true.
      bool weed_out(); // removes unconvex points one by one
  };
  //! Prints out convex-hull using usual << operator.
  template< class T_OBJECT >
  std::ostream& operator<<( std::ostream &_os, const Base<T_OBJECT> &_b )
  {
    return _os << _b.print();
  }
      
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
       Print :: out << "Storing first vertex " << vertex 
                    << *this << Print::endl;
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
       Print :: out << "Storing second vertex " << vertex 
                    << *this << Print::endl;
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
     
       Print :: out << "Storing other vertex " << vertex
                    << *this << Print::endl;
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

    std::sort( sortcopy.begin(), sortcopy.end(), std::mem_fun_ref( &Vertex<t_Object>::x_less ) );
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
  bool Base<T_OBJECT> :: Save( TiXmlElement &_node, T_SAVEOP &_op ) const
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
  bool Vertex<T_OBJECT> :: Save( TiXmlElement &_node, T_SAVEOP &_op ) const
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
}

#endif //  _CONVEX_HULL_BASE_H_
