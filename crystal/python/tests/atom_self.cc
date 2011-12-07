#include "LaDaConfig.h"

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/object.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/reference_existing_object.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/ref.hpp>

#include <python/exceptions.h>
#include "../atom/atom.h"

namespace bp = boost::python;
using namespace LaDa;
#if LADA_TYPE == 0
  typedef crystal::Atom<std::string> Atom;
# define NewAtom PyAtomStr_New
# define LADA_MODULE atom_self_0
#elif LADA_TYPE == 1 
  typedef crystal::Atom< std::vector<std::string> > Atom;
# define NewAtom PyAtomSequence_New
# define LADA_MODULE atom_self_1
#endif

bp::object get_new_object() { return bp::object(bp::handle<>(NewAtom())); }
bp::object get_static_object()
{
  static Atom atom;
  return bp::object(bp::handle<>(PyAtom_FromAtom(atom))); 
}

BOOST_PYTHON_MODULE(LADA_MODULE)
{
  bp::def("get_new_object", &get_new_object);
  bp::def("get_static_object", &get_static_object);
}
#undef NewAtom
#undef LADA_MODULE
