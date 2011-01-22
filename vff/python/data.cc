#include "LaDaConfig.h"

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/make_constructor.hpp>


#include <python/numpy_types.h>

#include "data.hpp"
#include "../data.h"

namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;
    namespace mn = math::numpy;
    vff::BondData* create_BondData()
    {
      try
      {
        vff::BondData *result = new vff::BondData;
        result->length = 0;
        for(size_t i(0); i < vff::max_vff_expansion; ++i) result->alphas[i] = 0e0;
        return result;
      }
      catch(std::exception &_e)
      {
        PyErr_SetString(PyExc_RuntimeError, "Could not create BondData object.");
        bp::throw_error_already_set();
      }
      return NULL;
    }

    vff::AngleData* create_AngleData()
    {
      try
      {
        vff::AngleData *result = new vff::AngleData;
        result->gamma = 0e0;
        result->sigma = 0e0;
        for(size_t i(0); i < vff::max_vff_expansion; ++i) result->betas[i] = 0e0;
        return result;
      }
      catch(std::exception &_e)
      {
        PyErr_SetString(PyExc_RuntimeError, "Could not create AngleData object.");
        bp::throw_error_already_set();
      }
      return NULL;
    }
    
    bp::object get_alphas(vff::BondData &_this)
    {
      std::cout << "THERE " << vff::max_vff_expansion << "\n";
      bp::object result = mn::create_array<mn::type<types::t_real>::value>(vff::max_vff_expansion); //_this.alphas, vff::max_vff_expansion);
//     std::cout << "HERE" << PyArray_REFCOUNT(result.ptr()) << "\n" << std::flush;
      return bp::object();
    }
    void set_alphas(vff::BondData &_this, bp::object _object)
    {
      std::ostringstream sstr; sstr << vff::max_vff_expansion; 
      size_t const N(bp::len(_object));
      if(N > vff::max_vff_expansion)
      {
        PyErr_SetString(PyExc_ValueError, ("Cannot set alphas beyond order " + sstr.str() + ".").c_str());
        bp::throw_error_already_set();
        return;
      }
      std::cout << "There " << _this.alphas << std::flush;
      for(size_t i(0); i < vff::max_vff_expansion; ++i)
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
    inline bp::object get_betas(vff::AngleData &_this)
    {
      return math::numpy::array_from_ptr(_this.betas, vff::max_vff_expansion);
    }
    inline void set_betas(vff::AngleData &_this, bp::object _object)
    {
      std::ostringstream sstr; sstr << vff::max_vff_expansion; 
      size_t const N(bp::len(_object));
      if(N > vff::max_vff_expansion)
      {
        PyErr_SetString(PyExc_ValueError, ("Cannot set betas beyond order " + sstr.str() + ".").c_str());
        bp::throw_error_already_set();
        return;
      }
      for(size_t i(0); i < vff::max_vff_expansion; ++i)
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

    void expose_data()
    {
      std::ostringstream sstr; sstr << vff::max_vff_expansion;
      bp::class_<vff::BondData>("BondData", "Bond-parameters for vff.")
        .def("__init__", bp::make_constructor(&create_BondData) )
        .def_readwrite("length", &vff::BondData::length, "Equilibrium bond-length")
        .add_property( "alphas", 
                       bp::make_function(&get_alphas, bp::with_custodian_and_ward_postcall<0, 1>()),
                       bp::make_function(&set_alphas));
      bp::class_<vff::AngleData>("AngleData", "Angle-parameters for vff.")
        .def("__init__", bp::make_constructor(&create_AngleData) )
        .def_readwrite("sigma", &vff::AngleData::sigma, "Equilibrium bond-length")
        .def_readwrite("gamma", &vff::AngleData::gamma, "Equilibrium bond-length")
        .add_property( "betas", 
                       bp::make_function(&get_betas, bp::with_custodian_and_ward_postcall<0, 1>()),
                       &set_betas);
      bp::def("_bond_type", &vff::bond_type);
      bp::def("_angle_type", &vff::bond_type);
    }

  }
} // namespace LaDa
