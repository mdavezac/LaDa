#include "LaDaConfig.h"

#include <vector>

#include <algorithm>
#include <numeric>
#include <functional>

#include <boost/lambda/lambda.hpp>
#include <boost/python/def.hpp> 
#include <boost/python/tuple.hpp> 
#include <boost/python/object.hpp> 

#include <opt/types.h>
#include <minimizer/python/minimizer.hpp>
#include "../functional.h"
#include "functional.hpp"

namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;

    bp::tuple minimize( Models::Clj const & _functional,
                        bp::object &_minimizer,
                        Models::Clj::t_Arg &_structure,
                        bool relax_all )
    {
      try
      {
        namespace bl = boost::lambda;
        Models::Functional functional( _functional );
        Models::Functional::t_Arg args;
        functional.relaxation = relax_all ? Models::Functional::Relaxation::default_:
                                            Models::Functional::Relaxation::volume;  
        functional.structure = _structure;
        functional.init( args );

        math::rVector3d const center 
        (
          std::accumulate
          (
             _structure.atoms.begin(), _structure.atoms.end(),
             math::rVector3d(0,0,0), std::plus<math::rVector3d>()
          )
        );

        bp::object object( _minimizer( functional, args ) ); 
        Models::Functional::t_Return result = bp::extract<Models::Functional::t_Return>(object);

        math::rVector3d const center2 
        (
          (  
            center - std::accumulate
                     ( 
                       functional.structure.atoms.begin(), functional.structure.atoms.end(), 
                       math::rVector3d(0,0,0), std::plus<math::rVector3d>()
                     )
          ) / types::t_real( _structure.atoms.size() )
        );

        foreach( Models::Clj::t_Arg::t_Atom& atom, functional.structure.atoms ) atom.pos -= center2;


        _structure.energy = result;
        functional.forces.energy = result;
        return bp::make_tuple(functional.forces, functional.structure);
      }
      catch (...)
      {
        PyErr_SetString( PyExc_RuntimeError, "Error encountered while minimizing from C.\n");
        return bp::tuple();
      }
    }
  
    void expose_functional()
    {
      bp::def
      ( 
        "minimize", 
        &minimize,
        (
          bp::arg("clj"),
          bp::arg("minimizer"),
          bp::arg("structure"),
          bp::arg("relaxation") = true
        ),
        "Minimizes a Coulomb +Lennard-Jhones functional.\n\n"
        ":Parameter:\n"
        "clj : `lada.pcm.Clj`\n  Coulomb + Lennard-Jones functional.\n"
        "minimizer : `lada.minimizer.Minimizer`\n   Minimization object.\n"
        "structure : `lada.crystal.Structure` \n  Input crystal structure.\n"
        "relaxation \n  If True, relax everything. If False only relax volume.\n"
        ":return: an 2-tuple of `lada.crystal.Structure` objects containing the "
        "relaxed structure for the first, and the forces and stress for the "
        "second."
      );
    }
  } // namespace python
} // namespace LaDa
