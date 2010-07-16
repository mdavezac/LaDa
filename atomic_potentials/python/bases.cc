#include "LaDaConfig.h"

#include <sstream>
#include <complex>

#include <boost/python/class.hpp>
#include <boost/python/errors.hpp>

#include "../bases.h"


namespace LaDa
{
  namespace Python
  {

    struct Bases 
    {
      typedef atomic_potential::Bases<Crystal::TStructure<std::string> >::const_iterator t_cit;
      Bases   ( Crystal::TStructure<std::string> const &_str )
            : bases_(_str), first_(true)
        { cit_ = bases_.begin(); cit_end_ = bases_.end(); }
      Bases   ( Bases const &_c )
            : bases_(_c.bases_), cit_(_c.cit_), cit_end_(_c.cit_end_), first_(_c.first_) {}

      Bases &iter()  { return *this; }
      t_cit::reference next()
      {
        namespace bp = boost::python;
        if( first_ ) first_ = false; 
        else 
        {
          t_cit const copy(cit_);
          ++cit_;
          if( cit_ == cit_end_ )
          {
            PyErr_SetString
            (
              PyExc_StopIteration, 
              "Error while computing transform to smith normal form.\n" 
            );
            bp::throw_error_already_set();
           cit_ = copy;
          }
        }
        return *cit_;
      }

      atomic_potential::Bases<Crystal::TStructure<std::string> > bases_;
      t_cit cit_;
      t_cit cit_end_;
      bool first_;
    };

    void expose_bases()
    {
      namespace bp = boost::python;

      typedef atomic_potential::Basis Basis;
      bp::class_<Basis>("Basis", "Cartesian basis.")
          .add_property("origin", &Basis::origin, "Origin of the basis.")
          .add_property("x", &Basis::x, "x unit vector.")
          .add_property("y", &Basis::y, "y unit vector.")
          .add_property("z", &Basis::z, "z unit vector.")
          .add_property("weight", &Basis::weight, "weight of the basis in a representation.")
          .add_property("index", &Basis::weight, "Index of the origin in the structure.");

      bp::class_<Bases>
      (
        "Bases", 
        "Iterator over the basis of a representation.\n",
        bp::init<Crystal::TStructure<std::string> const &>()
      ).def("__iter__", &Bases::iter, bp::return_internal_reference<1>())
       .def("next", &Bases::next, bp::return_internal_reference<1>());
    }

  }
} // namespace LaDa
//  Version: $Id$
//

#ifndef LADA_CRYSTAL_PYTHON_NEIGHBORS_HPP
#define LADA_CRYSTAL_PYTHON_NEIGHBORS_HPP

#include "LaDaConfig.h"

namespace LaDa
{
  namespace Python
  {
    void expose_neighbors();
  }
} // namespace LaDa
#endif
