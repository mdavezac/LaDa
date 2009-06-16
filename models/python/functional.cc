//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <algorithm>
#include <numeric>
#include <functional>

#include <boost/lambda/lambda.hpp>
#include <boost/python/def.hpp> 
#include <boost/python/tuple.hpp> 
#include <boost/python/enum.hpp> 

#include <opt/types.h>
#include <minimizer/python/minimizer.hpp>
#include "../functional.h"
#include "functional.hpp"

namespace LaDa
{
  namespace Python
  {
    Models::Clj::t_Arg minimize( Models::Clj const & _functional,
                                 Python::Minimizer const &_minimizer,
                                 Models::Clj::t_Arg &_structure,
                                 Models::Functional::Relaxation::type _relaxation )
    {
      try
      {
        namespace bl = boost::lambda;
        Models::Functional functional( _functional );
        Models::Functional::t_Arg args;
        functional.relaxation = _relaxation;
        functional.structure = _structure;
        functional.init( args );

        atat::rVector3d const center 
        (
          std::accumulate
          (
             _structure.atoms.begin(), _structure.atoms.end(),
             atat::rVector3d(0,0,0), std::plus<atat::rVector3d>()
          )
        );

        Models::Functional::t_Return result( _minimizer.call( functional, args ) ); 

        atat::rVector3d const center2 
        (
          (  
            center - std::accumulate
                     ( 
                       functional.structure.atoms.begin(), functional.structure.atoms.end(), 
                       atat::rVector3d(0,0,0), std::plus<atat::rVector3d>()
                     )
          ) / types::t_real( _structure.atoms.size() )
        );

        foreach( Models::Clj::t_Arg::t_Atom& atom, functional.structure.atoms ) atom.pos -= center2;


        _structure.energy = result;
        functional.forces.energy = result;
        _structure = functional.structure;
        return functional.forces;
      }
      catch (...)
      {
        PyErr_SetString( PyExc_RuntimeError, "Error encountered while minimizing from C.\n");
        Models::Clj::t_Arg forces;
        forces.cell.zero();
        forces.atoms.clear();
        forces.energy = 0e0;
        return forces;
      }
    }
  
    void expose_functional()
    {
      namespace bp = boost :: python;
      bp::enum_< Models::Functional::Relaxation::type >( "relaxation" )
        .value( "default", Models::Functional::Relaxation::default_ )
        .value( "volume",  Models::Functional::Relaxation::volume )
        .export_values();

      bp::def
      ( 
        "minimize", 
        &minimize,
        (
          bp::arg("clj"),
          bp::arg("minimizer"),
          bp::arg("structure"),
          bp::arg("relaxation") = Models::Functional::Relaxation::default_
        ),
        "Minimizes a Coulomb +Lennard-Jhones functional:\n"
        "  clj is the functional.\n"
        "  minimizer is the minimizer object.\n"
        "  structure is the input/output cell and atomic-position.\n"
        "  relaxation can be either models.relaxation.default, or models.relaxation.volume.\n"
        "The return is a structure object containing the stress and the forces."
      );
    }
  } // namespace Python
} // namespace LaDa
