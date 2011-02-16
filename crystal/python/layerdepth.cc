#include "LaDaConfig.h"

#include <sstream>
#include <complex>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/init.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

#include <math/fuzzy.h>
#include "../atom.h"
#include "../structure.h"
#include "../layerdepth.h"

#include "layerdepth.hpp"


namespace LaDa
{
  namespace Python
  {
    namespace bp = boost::python;
    template<class T> inline math::rVector3d const& get_pos(T const &_a) { return _a.pos; }
    template<> inline math::rVector3d const& get_pos( math::rVector3d const &_a) { return _a; }

    template<class T1, class T2>
      bp::object call( Crystal::LayerDepth const& _l, T1 const &_a, T2 const &_b )
      {
        bool const sup( _l( get_pos(_a), get_pos(_b) ) );
        bool const inf( _l( get_pos(_b), get_pos(_a) ) );
        if( (not sup) and (not inf) ) return bp::object(0);
        return bp::object(sup ? -1: 1);
      }
    template<class T>  types::t_real call_single( Crystal::LayerDepth const &_l, T const &_a )
      { return _l(get_pos(_a)); }

    template<class T_TYPE>
      boost::shared_ptr< T_TYPE > sort_structure2( T_TYPE &_structure, math::rVector3d const &_dir )
      {
        boost::shared_ptr< T_TYPE > result(new T_TYPE(_structure));
        std::sort( result->atoms.begin(), result->atoms.end(), Crystal::LayerDepth(_dir) );
        return result;
      }
    template<class T_TYPE>
      boost::shared_ptr< T_TYPE > sort_structure( T_TYPE &_structure )
        { return sort_structure2(_structure, _structure.cell.col(0)); }

    void expose_layerdepth()
    {
      typedef math::rVector3d t_vec;
      typedef bp::return_value_policy<bp::return_by_value> t_repol;
      bp::class_<Crystal::LayerDepth>
      ( 
        "LayerDepth", 
        "Strict Weak Ordering functor according to depth along eptiaxial direction\n\n"
        "Computes the superlattice z-order of a given vector or compares the "
        "z-order of two vectors.",
        bp::no_init 
      ).def( bp::init<math::rVector3d const&>() )
       .def( bp::init<Crystal::LayerDepth const&>() )
       .def( "__call__", &call_single<t_vec> )
       .def
       (
         "__call__", 
         &call<t_vec, t_vec>,
         "If one argument, returns layer-depth. If two arguments, returns comparison.\n\n"
         "Arguments are numpy 3x1 float arrays."
       );
      bp::def( "sort_layers", &sort_structure<Crystal::Structure>, bp::arg("structure") );
      bp::def( "sort_layers", &sort_structure< Crystal::TStructure<std::string> >,
                bp::arg("structure") );
      bp::def( "sort_layers", &sort_structure2<Crystal::Structure>,
               (bp::arg("structure"), bp::arg("direction")) );
      bp::def
      ( 
        "sort_layers", 
        &sort_structure2< Crystal::TStructure<std::string> >,
        (bp::arg("structure"), bp::arg("direction")), 
        "Returns a structure with the atoms sorted by layer.\n\n"
        ":Parameters:\n"
        "  structure\n    An input structure. Not modified on output.\n"
        "  direction\n    Optional argument giving the direction and length of "
           "repetitions. When not given, this is the firstr cell vector of the structure, "
           "``structure.cell[:,0]``.\n"
      );
    }

  }
} // namespace LaDa
