#include "LaDaConfig.h"

#include <sstream>
#include <utility>

#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/def.hpp>
#include <boost/python/errors.hpp>
#include <boost/shared_ptr.hpp>

#define LADA_PYTHON_STD_VECTOR_NOPRINT
#include <python/std_vector.hpp>

#include "periodic_dnc.hpp"
#include "../periodic_dnc.h"


namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    typedef Crystal::DnCBoxes DnCBoxes;
    typedef boost::shared_ptr<DnCBoxes> ptr_DnCBoxes;
    typedef Crystal::DnCBoxes::t_Box DnCBox;
    typedef std::pair<DnCBoxes::t_Boxes::const_iterator,
                      DnCBoxes::t_Boxes::const_iterator> t_BoxesPair;

    struct PointIterator 
    {
      typedef Crystal::DnCBoxes::t_Box DnCBox;
      typedef DnCBox::const_iterator const_iterator;
      ptr_DnCBoxes ptr_boxes;
      std::pair<const_iterator, const_iterator> interval;
      bool is_first;
      PointIterator () : ptr_boxes() {}
      PointIterator   (ptr_DnCBoxes const &_boxes, t_BoxesPair::first_type const &_first)
                    : ptr_boxes(_boxes), interval(_first->begin(), _first->end()),
                      is_first(true) {}
      PointIterator   (PointIterator const & _c)
                    : ptr_boxes(_c.ptr_boxes), interval(_c.interval), is_first(_c.is_first) {};
      bp::tuple next() 
      {
        if(not ptr_boxes)
        {
          PyErr_SetString(PyExc_StopIteration, "PointIterator");
          bp::throw_error_already_set();
          return bp::tuple();
        }
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
    };

    struct BoxIterator 
    {
      ptr_DnCBoxes ptr_boxes;
      t_BoxesPair interval;
      bool is_first;
      BoxIterator( Crystal::TStructure<std::string> const &_structure,
                   int _size, types::t_real _overlap ) : is_first(true)
      {
        try
        { 
          ptr_boxes = ptr_DnCBoxes(new Crystal::DnCBoxes);
          math::iVector3d const mesh( ptr_boxes->guess_mesh(_structure, _size) );
          ptr_boxes->init(_structure, mesh, _overlap);
          interval.first = ptr_boxes->begin(); 
          interval.second = ptr_boxes->end();
        } catch(std::exception &_e)
        {
          std::ostringstream sstr;
          sstr << "Encountered error while creating DnCBoxes.\n"
               << _e.what() << "\n";
          PyErr_SetString(PyExc_RuntimeError, sstr.str().c_str());
          bp::throw_error_already_set();
        }
      }
      BoxIterator   (BoxIterator const &_c)
                  : ptr_boxes(_c.ptr_boxes), interval(_c.interval), is_first(_c.is_first) {}
      PointIterator next() 
      {
        if(is_first) is_first = false;
        else if(interval.first != interval.second) ++interval.first;
        if(interval.first == interval.second) 
        {
          PyErr_SetString(PyExc_StopIteration, "");
          bp::throw_error_already_set();
          return PointIterator();
        }
        return PointIterator(ptr_boxes, interval.first);
      }
      math::iVector3d mesh() const { return ptr_boxes->mesh(); }
      types::t_real overlap() const { return ptr_boxes->overlap(); }
    };

    inline bp::object pass_through(bp::object const& o) { return o; }


    void expose_periodic_dnc()
    {
      bp::class_<BoxIterator>
        (
          "dnc_iterator", 
          "Divide and Conquer mesh of boxes for periodic systems.\n\n"
          ":Parameters:\n"
          "  structure : `Structure`\n    Periodic structure to divide in real-space.\n"
          "  size : int\n    Tentative number of atoms per box.\n"
          "  overlap : float\n    Overlap between boxes, in angstrom."
          "A DnCBoxes mesh will parcel out all atoms in a structure, "
          "including periodic images. In addition, boxes are defined with an overlap, "
          "so that atom at the outside edge of box are also included in the box. "
          "In pratice, should be used as an iterator over the atoms. ",
          bp::no_init
        )
        .def(bp::init< Crystal::TStructure<std::string> const&, int, types::t_real>() )
        .def_readonly("mesh", &BoxIterator::mesh)
        .def_readonly("overlap", &BoxIterator::overlap)
        .def("__iter__", &pass_through)
        .def("next", &BoxIterator::next);

      bp::class_<PointIterator>("DnCBoxIterator", bp::no_init)
        .def("__iter__", &pass_through)
        .def("next", &PointIterator::next);
    }

  } // namespace Python
} // namespace LaDa
