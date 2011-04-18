#include "LaDaConfig.h"

#include <sstream>
#include <complex>

#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/def.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/return_internal_reference.hpp>

#define LADA_PYTHON_STD_VECTOR_NOPRINT
#include <python/std_vector.hpp>

#include "periodic_dnc.hpp"
#include "../periodic_dnc.h"


namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    struct Iterator 
    {
      std::pair<Crystal::DnCBox::const_iterator, Crystal::DnCBox::const_iterator> interval;
      bool is_first;
      Iterator(Crystal::DnCBox const &_box) : interval(_box.begin(), _box.end()), is_first(true) {}
      Crystal::DnCBox::const_iterator::value_type next() 
      {
        if(is_first) is_first = false;
        else if(interval.first != interval.second) ++interval.first;
        if(interval.first == interval.second) 
        {
          PyErr_SetString(PyExc_StopIteration, "");
          bp::throw_error_already_set();
          static Crystal::DnCBox::const_iterator::value_type dummy;
          return dummy;
        }
        return bp::make_tuple(interval.first->index, interval.first->trans, interval.first->end());
      }
      Iterator & iter() const { return *this; }
    }
    Iterator create(Crystal::DnCBox const &_box) { return Iterator(_box); }

    void expose_periodic_dnc()
    {
      bp::def( "create_periodic_dnc", &Crystal::create_periodic_dnc<std::string>,
               (bp::arg("structure"), bp::arg("overlap")), "Creates a list of dnc boxes." );
      bp::class_<Crystal::DnCBox>
        ( "DnCBox", bp::noinit,
          "Divide and Conquer box for periodic systems.\n\n"
          "A mesh DnCBox objects will parcel out all atoms in a structure, "
          "including periodic images. In addition, boxes are defined with an overlap, "
          "so that atom at the outside edge of box are also included in the box. "
          "In pratice, should be used as an iterator over the atoms. " )
        .def( "__iter__", &create,
              "Iterates over atoms in/around DnCBox.\n\n"
              "Yields a 3-tuple consisting of the index to the atom in the structure, "
              "a vector acounting for possible periodic translations, "
              "a boolean which is true if the atom is in the box (as opposed to its outside edge." );
      bp::class_<Iterator>("DnCBoxIterator", bp::noinit)
        .def("__iter__", &Iterator::iter, bp::return_internal_reference<1>())
        .def("next", &Iterator::next);
      expose_vector<Crystal::t_DnCBoxes::type>("DnCBoxes", "Mesh of divide and conquer boxes.");
      bp::register_ptr_to_python<Crystal::t_DnCBoxes::shared_ptr>();
    }

  } // namespace Python
} // namespace LaDa
