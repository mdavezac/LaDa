#include "LaDaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/object.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/reference_existing_object.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/ref.hpp>

#include <python/exceptions.h>
#include "../../atom.h"
#include "../../traits.h"

namespace bp = boost::python;
using namespace LaDa;
typedef crystal::traits::StructureData< LADA_TYPE >::t_Atom Atom;

bp::object get_new_object() { return bp::object(Atom()); }
Atom get_new_atom() { return Atom(); }
Atom& get_incorrect_reference()
{
  static Atom atom;
  return atom;
}
bp::object get_reference()
{
  static Atom atom;
  if(atom.self().ptr() != Py_None) return atom.self();
  bp::reference_existing_object::apply<Atom const*>::type converter;
  bp::object o = bp::object(bp::handle<>(converter(&atom)));
  o.attr("_set_self")();
  return o;
}
bp::object get_pointer()
{
  std::auto_ptr<Atom> atom(new Atom);
  bp::object o = bp::object(atom);
  o.attr("_set_self")();
  return o;
}

void check_back_reference(bp::back_reference<Atom const&> _ref, Atom const &_ob)
{
  if(_ref.get().self() != _ref.source())
    python::PyException<error::InternalError>::throw_error("back-ref's source and self_ not equivalent.");
  if(_ref.get().self() != _ref.source())
    python::PyException<error::InternalError>::throw_error("back-ref's source and atom's self_ not equivalent.");
}
void change_pos(Atom &_ob) { _ob.pos[0] = 1; }


BOOST_PYTHON_MODULE(LADA_MODULE)
{
  bp::def("get_new_object", &get_new_object);
  bp::def("get_new_atom", &get_new_atom);
  bp::def("get_reference", &get_reference);
  bp::def("get_pointer", &get_pointer);
  bp::def( "get_incorrect_reference", &get_incorrect_reference,
           bp::return_value_policy<bp::reference_existing_object>());
  bp::def( "check_back_reference", &check_back_reference);
  bp::def( "change_pos", &change_pos);
}
