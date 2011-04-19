#include "LaDaConfig.h"

#include <sstream>
#include <utility>

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
    struct BoxIterator 
    {
      typedef Crystal::DnCBoxes::t_Boxes::const_iterator const_iterator;
      std::pair<const_iterator, const_iterator> interval;
      bp::object object;
      bool is_first;
      BoxIterator   (Crystal::DnCBoxes const &_boxes, bp::object const &_object)
                  : interval(_boxes.begin(), _boxes.end()), object(_object), is_first(true) {}
      BoxIterator   (Crystal::DnCBoxes const &_boxes, bp::object const &_object)
                  : interval(_boxes.begin(), _boxes.end()), object(_object), is_first(true) {}
      const_iterator::value_type const & next() 
      {
        if(is_first) is_first = false;
        else if(interval.first != interval.second) ++interval.first;
        if(interval.first == interval.second) 
        {
          PyErr_SetString(PyExc_StopIteration, "");
          bp::throw_error_already_set();
          const_iterator::value_type const static dummy;
          return dummy; 
        }
        return *interval.first;
      }
      BoxIterator & iter() { return *this; }
    };
    BoxIterator create_boxiterator(Crystal::DnCBoxes const &_box)
      { return BoxIterator(_box, bp::object()); }

    BoxIterator dnx_iterator( Crystal::TStructure<std::string> const &_structure,
                              int _size, types::t_real _overlap)
    {
      bp::object result( bp::handle<>(Crystal::DnCBoxes()) );
      Crystal::DnCBoxes &boxes = bp::extract<Crystal::DnCBoxes &>(result);
      math::iVector3d const mesh( boxes.guess_mesh(structure) );
      boxes.init(mesh, _overlap);
      return BoxIterator(boxes);
    }

    struct PointIterator 
    {
      typedef Crystal::DnCBoxes::t_Box DnCBox;
      typedef DnCBox::const_iterator const_iterator;
      std::pair<const_iterator, const_iterator> interval;
      bool is_first;
      PointIterator(DnCBox const &_box) : interval(_box.begin(), _box.end()), is_first(true) {}
      bp::tuple next() 
      {
        if(is_first) is_first = false;
        else if(interval.first != interval.second) ++interval.first;
        if(interval.first == interval.second) 
        {
          PyErr_SetString(PyExc_StopIteration, "");
          bp::throw_error_already_set();
          return bp::tuple();
        }
        return bp::make_tuple( interval.first->index,  
                               interval.first->translation, 
                               interval.first->in_small_box );
      }
      PointIterator & iter() { return *this; }
    };
    PointIterator create(Crystal::DnCBoxes::t_Box const &_box) { return PointIterator(_box); }

    inline bp::object pass_through(bp::object const& o) { return o; }


    void expose_periodic_dnc()
    {
      typedef bp::return_value_policy< bp::with_custodian_and_ward_postcall<1, 0> > t_rvp;
      bp::def( "dnc_iterator", &dnc_iterator,
               (bp::arg("structure"), bp::arg("size"), bp::arg("overlap")=1.5),
               "Iterates over Divide-and-Conquer boxes.\n\n"
               ":Parameters:\n"
               "  structure : `Structure`\n    Periodic structure to divide in real-space.\n"
               "  size : int\n    Tentative number of atoms per box.\n"
               "  overlap : float\n    Overlap between boxes, in angstrom." );
                 

      bp::class_<Crystal::DnCBoxes>
        ( "DnCBoxes", 
          "Divide and Conquer mesh of boxes for periodic systems.\n\n"
          "A DnCBoxes mesh will parcel out all atoms in a structure, "
          "including periodic images. In addition, boxes are defined with an overlap, "
          "so that atom at the outside edge of box are also included in the box. "
          "In pratice, should be used as an iterator over the atoms. " )
        .def_readonly( "overlap",   &Crystal::DnCBoxes::overlap )
        .def_readonly( "mesh",   &Crystal::DnCBoxes::mesh )
        .def("__iter__", &create_boxiterator);

      bp::class_<Crystal::DnCBoxes::t_Box>
        ( "DnCBox", bp::no_init )
        .def( "__iter__", &create,
              "Iterates over atoms in/around a box.\n\n"
              "Yields a 3-tuple consisting of the index to the atom in the structure, "
              "a vector acounting for possible periodic translations, "
              "a boolean which is true if the atom is in the box (as opposed to its outside edge." )
        .def("__iter__", &create);
      bp::class_<PointIterator>("DnCBoxIterator", bp::no_init)
        .def("__iter__", &pass_through)
        .def("next", &PointIterator::next);
    }

  } // namespace Python
} // namespace LaDa
