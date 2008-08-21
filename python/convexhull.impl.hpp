//
//  Version: $Id$
//

#ifndef _INMODULE_


#ifndef _TYPE_
#  error "Define macro _TYPE_ to use this file."
#endif
#ifndef _NAMESPACE_
#  error "Define macro _NAMESPACE_ to use this file."
#endif
namespace _NAMESPACE_
{
  typedef ::opt::ConvexHull< ::opt::ConvexHull::FakeObject< _TYPE_ > > t_FakeObject;
  typedef ::opt::ConvexHull::Vertex< t_FakeObject > t_VertexObject;
  typedef ::opt::ConvexHull::Base< t_VertexObject > t_CHBase;

  void add_to_ch( t_CHBase &_ch, const tuple &_t )
  {
    __DOASSERT( len( _t ) != 3, "Convex-hull tuple must contain three objects." )
    t_FakeObject fake( extract< _TYPE_ >( _t[0] ) );
   ::opt::ConvexHull::NullObject::x = extract< types::t_real >( _t[1] );
    types::t_real energy = extract< types::t_real >( _t[2] );
    _ch.add( fake, energy );
  }
}
#undef _TYPE_
#undef _NAMESPACE_


#else 


#ifndef _PYTHONNAME_
#  error "Define macro _PYTHONNAME_ to use this file."
#endif
#ifndef _NAMESPACE_
#  error "Define macro _NAMESPACE_ to use this file."
#endif

  // Exposes convex-hull.
  class_< _NAMESPACE_::t_CHBase >( _PYTHONNAME_ )
    .def( "evaluate", &_NAMESPACE_::t_CHBase::evaluate )
    .def( "add", &_NAMESPACE_::add_to_ch )
    .def( "__str__", &_NAMESPACE_::t_CHBase::print );

#undef _INMODULE_
#undef _PYTHONNAME_
#undef _NAMESPACE_



#endif
