#include "LaDaConfig.h"

#include <sstream>
#include <complex>

#include <boost/python/class.hpp>
#include <boost/python/errors.hpp>

#include "../neighbors.h"


namespace LaDa
{
  namespace Python
  {

    struct Neighbors 
    {
      Neighbors   ( Crystal::TStructure<std::string> const &_str,
                     size_t _nmax = 0, 
                     math::rVector3d const &_orig = math::rVector3d(0,0,0) )
                 : neighs_(_nmax, _orig), first_(true)
        { cit_ = neighs_.begin(_str); cit_end_ = neighs_.end(); }
      Neighbors   ( Neighbors const &_c )
                 : neighs_(_c.neighs_), cit_(_c.cit_), cit_end_(_c.cit_end_), first_(_c.first_) {}

      Neighbors &iter()  { return *this; }
      Crystal::Neighbors::const_iterator::reference next()
      {
        namespace bp = boost::python;
        if( first_ ) first_ = false; 
        else 
        {
          ++cit_;
          if( cit_ == cit_end_ )
          {
            PyErr_SetString
            (
              PyExc_StopIteration, 
              "Error while computing transform to smith normal form.\n" 
            );
            bp::throw_error_already_set();
            --cit_;
          }
        }
        return *cit_;
      }

      Crystal::Neighbors neighs_;
      Crystal::Neighbors::const_iterator cit_;
      Crystal::Neighbors::const_iterator cit_end_;
      bool first_;
    };

    void expose_neighbors()
    {
      namespace bp = boost::python;

      typedef Crystal::Neighbor Neighbor;
      bp::class_<Neighbor>("Neighbor", "Describes nth neighbor.\n\n"
                           "This is the object returned by L{Neighbors}.")
          .def(bp::init<Neighbor const&>())
          .add_property("index", &Neighbor::index, "Index of atom in structure.")
          .add_property("pos", &Neighbor::pos, "Position respect to the origin.")
          .add_property("distance", &Neighbor::distance, "Distance from the origin.");

      bp::class_<Neighbors>
      (
        "Neighbors", 
        "Iterator over first neighbors.\n\n"
        "Initialization:\n"
        " - The first argument (L{Structure}) is the structure within which to find neighbors.\n"
        " - The second argument (integer) is the number of neighbors to look for.\n"
        " - The third argument (3x numpy array) is the central atom position "
        "for which to find neightbors.",
        bp::init<Crystal::TStructure<std::string> const &, size_t, math::rVector3d>()
      ).def("__iter__", &Neighbors::iter, bp::return_internal_reference<1>())
       .def("next", &Neighbors::next, bp::return_internal_reference<1>());
    }

  }
} // namespace LaDa
