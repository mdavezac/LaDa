//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/python/class.hpp>
#include <boost/python/scope.hpp>

#include "values.h"
#include "../collapse/values.h"
#include "../sum_of_separables.h"
#include "../representation.h"
#include "../collapse/fitting_set.h"
#define LADA_PYTHON_STD_VECTOR_NOPRINT
#include <python/std_vector.hpp>


namespace LaDa
{
  namespace Python
  {
    typedef LaDa::atomic_potential::collapse::Values Values;
    atomic_potential::numeric_type other( Values const& _val, size_t _i, size_t _r,
                                          size_t _s, size_t _rep )
    {
      namespace bp = boost::python;
      if( _i >= _val.nb_coordinates() )
      {
        PyErr_SetString(PyExc_IndexError, "Coordinate index out of range.\n");
        bp::throw_error_already_set();
        return -1;
      }
      Values::const_str_iterator i_str = _val.begin(_i);
      Values::const_str_iterator const i_str_end = _val.end(_i);
      for(size_t s(_s);  i_str != i_str_end and s != 0; --s, ++i_str);
      if( i_str == i_str_end )
      {
        PyErr_SetString(PyExc_IndexError, "Structure index out of range.\n");
        bp::throw_error_already_set();
        return -1;
      }
      Values::const_str_iterator::const_rep_iterator i_rep = i_str.begin();
      Values::const_str_iterator::const_rep_iterator const i_rep_end = i_str.end();
      for(size_t rep(_rep);  i_rep != i_rep_end and rep != 0; --rep, ++i_rep);
      if( i_rep == i_rep_end )
      {
        PyErr_SetString(PyExc_IndexError, "Representation index out of range.\n");
        bp::throw_error_already_set();
        return -1;
      }
      typedef Values::const_str_iterator::const_rep_iterator::const_rank_iterator t_cit;
      t_cit i_rank = i_rep.begin();
      t_cit const i_rank_end = i_rep.end();
      for(size_t r(_r);  i_rank != i_rank_end and r != 0; --r, ++i_rank);
      if( i_rank == i_rank_end )
      {
        PyErr_SetString(PyExc_IndexError, "Rank index out of range.\n");
        bp::throw_error_already_set();
        return -1;
      }
      return i_rank.other();
    }
    void expose_values()
    {
      namespace bp = boost::python;
      
      bp::scope scope = bp::class_<Values>
            ("Values", "Holds fitting related info. For Debugging only.")
        .add_property("coord_rank_values", &Values::coord_rank_values_)
        .add_property("rank_values", &Values::rank_values_)
        .add_property("function_values", &Values::function_values_)
        .def
        (
          "other", 
          &other, 
          (
            bp::arg("coordinate"),
            bp::arg("rank"),
            bp::arg("structure"),
            bp::arg("representation")
          )
        )
        .def("add", &Values::add);

      expose_vector<Values::t_FunctionValues::value_type>
         ("Values_vec1", "Implementation only.");
      expose_vector<Values::t_FunctionValues::value_type::value_type>
         ("Values_vec2", "Implementation only.");
      expose_vector<Values::t_FunctionValues::value_type::value_type::value_type>
         ("Values_vec3", "Implementation only.");
      expose_vector<Values::t_FunctionValues::value_type::value_type::value_type::value_type>
         ("Values_vec4", "Implementation only.");
 
    }

  }
} // namespace LaDa
