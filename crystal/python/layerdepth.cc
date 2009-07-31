//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <sstream>
#include <complex>

#include <boost/python.hpp>
#include <boost/python/enum.hpp>

#include <opt/fuzzy.h>
#include "../atom.h"
#include "../structure.h"
#include "../layerdepth.h"

#include "layerdepth.hpp"


namespace LaDa
{
  namespace Python
  {
    types::t_int call( Crystal::LayerDepth const& _l,
                       atat::rVector3d const &_a,
                       atat::rVector3d const &_b )
    {
      bool const sup( _l(_a,_b) );
      bool const inf( _l(_b,_a) );
      if( (not sup) and (not inf) ) return 0;
      return sup ? -1: 1;
    }
    template<class T>
      types::t_int call_atom( Crystal::LayerDepth const& _l,
                              Crystal::Atom_Type<T> const &_a, 
                              Crystal::Atom_Type<T> const &_b )
        { return call( _l, _a.pos, _b.pos ); }

    template<class T_TYPE>
      void sort_structure( Crystal::TStructure<T_TYPE> &_structure )
      {
        std::sort( _structure.atoms.begin(), _structure.atoms.end(), 
                   Crystal::LayerDepth(_structure.cell) );
      }

    void expose_layerdepth()
    {
      namespace bp = boost::python;
      bp::class_<Crystal::LayerDepth>( "LayerDepth", bp::no_init )
        .def( bp::init<atat::rMatrix3d const&>() )
        .def( bp::init<Crystal::LayerDepth const&>() )
        .def( "__call__", &call_atom<std::string> )
        .def( "__call__", &call_atom<types::t_real> )
        .def( "__call__", &call );
      bp::def( "sort_layers", &sort_structure<std::string> );
      bp::def( "sort_layers", &sort_structure<types::t_real> );
    }

  }
} // namespace LaDa
