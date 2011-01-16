#include "LaDaConfig.h"

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>


#include <python/numpy_types.h>
#include <python/std_map.h>

#include "data.hpp"
#include "../data.h"

namespace LaDa
{
  namespace bp = boost::python;
  namespace python
  {
    vff::BondData create_BondData();
    {
      vff::BondData result;
      result.bond_length = 0;
      for(size_t i(0); i < vff::max_vff_expansion; ++i) result.alphas[i] = 0e0;
      return result;
    }

    vff::BondData create_AngleData();
    {
      vff::AngleData result;
      result.gamma = 0e0;
      result.sigma = 0e0;
      for(size_t i(0); i < vff::max_vff_expansion; ++i) result.betas[i] = 0e0;
      return result;
    }
    
    inline bp::object get_alphas(vff::BondData const &_this)
    {
      return math::numpy::array_from_ptr(_this.alphas, vff::max_vff_expansion);
    }
    inline void set_alphas(vff::BondData &_this, bp::object _object)
    {
      std::ostringstream sstr; sstr << vff::max_vff_expansion; 
      size_t const N(bp::len(_object));
      if(N > vff::max_vff_expansion)
      {
        PyErr_SetString(PyExc_ValueError, ("Cannot set alphas beyond order " + sstr.str() + ".").c_str());
        bp::throw_error_already_set();
        return;
      }
      for(size_t i(0); i < vff::max_vff_expension; ++i)
        if( i >= N ) _this.alphas[i] = 0e0;
        else
          try { _this.alphas[i] = bp::extract<types::t_real>(_object[i]); }
          catch(...)
          {
            PyErr_SetString(PyExc_ValueError, "Input value could not be converted to real number.");
            bp::throw_error_already_set();
            return;
          }
        }
    }
    inline bp::object get_betas(vff::BondData const &_this)
    {
      return math::numpy::array_from_ptr(_this.betas, vff::max_vff_expansion);
    }
    inline void set_betas(vff::BondData &_this, bp::object _object)
    {
      std::ostringstream sstr; sstr << vff::max_vff_expansion; 
      size_t const N(bp::len(_object));
      if(N > vff::max_vff_expansion)
      {
        PyErr_SetString(PyExc_ValueError, ("Cannot set betas beyond order " + sstr.str() + ".").c_str());
        bp::throw_error_already_set();
        return;
      }
      for(size_t i(0); i < vff::max_vff_expension; ++i)
        if( i >= N ) _this.betas[i] = 0e0;
        else
          try { _this.betas[i] = bp::extract<types::t_real>(_object[i]); }
          catch(...)
          {
            PyErr_SetString(PyExc_ValueError, "Input value could not be converted to real number.");
            bp::throw_error_already_set();
            return;
          }
        }
    }

    void expose_data()
    {
      std::ostringstream sstr; sstr << vff::max_vff_expansion;
      bp::class_<vff::BondData>("BondData", "Bond-parameters for vff.")
        .def("__init__", bp::make_constructor(&create_BondData) )
        .def_readwrite("length", "Equilibrium bond-length", bp::no_init)
        .add_property( "alphas", 
                       bp::make_function(&get_alphas, bp::with_custodian_and_ward<0, 1>()),
                       &set_alphas,
                       "Bond-stretching parameters.\n\n"
                       "First value is order two. Only up to order " + sstr.str() + " allowed." );
      bp::class_<vff::BondData>("BondData", "Bond-parameters for vff.")
        .def("__init__", bp::make_constructor(&create_AngleData) )
        .def_readwrite("length", "Equilibrium bond-length", bp::no_init)
        .add_property( "betas", 
                       bp::make_function(&get_betas, bp::with_custodian_and_ward<0, 1>()),
                       &set_betas,
                       "Bond-stretching parameters.\n\n"
                       "First value is order two. Only up to order " + sstr.str() + " allowed." );
      bp::def("_bond_type", &vff::bond_type);
      bp::def("_angle_type", &vff::angle_type);

      expose_map< std::map<std::string, vff::BondData> >("_BondDataMap");
      expose_map< std::map<std::string, vff::AngleData> >("_AngleDataMap");
    
    }

  }
} // namespace LaDa
