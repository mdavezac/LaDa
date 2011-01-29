#include "LaDaConfig.h"

#include <string>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/data_members.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/exception/get_error_info.hpp>
#include <boost/exception/diagnostic_information.hpp>
#ifdef LADA_MPI
# include <boost/mpi/collectives.hpp>
#endif


#include "vff.hpp"
#include "../functional.h"
#include "../layered.h"
#include "../with_minimizer.h"

namespace LaDa
{
  namespace bp = boost::python;
  namespace python
  {

    template<class T>
      void print_escan_input( T &_self, const std::string &_path, bp::object &_str )
      {
        if(_str.ptr() != Py_None)  
        {
          Crystal::TStructure<std::string> structure;
          try { structure = bp::extract< Crystal::TStructure<std::string> >(_str); }
          catch(std::exception &e)
          {
            std::ostringstream sstr;
            sstr << "Could not convert to structure.\n" << e.what() << "\n";
            PyErr_SetString(PyExc_ValueError, sstr.str().c_str());
            bp::throw_error_already_set();
            return;
          }
          _self.set_structure(structure);
          _self.init(true, false);
        }
        _self.print_escan_input(_path);
      }
    template<class T> 
      void init( T &_self, typename Crystal::TStructure<std::string> const &_str,
                 bool _doinit, bool _verbose = false )
      { 
        _self.VffBase().structure = _str;
        try { _self.init(_doinit, _verbose); }
        catch(vff::exceptions::site_index &e)
        {
          if(int const * mi = boost::get_error_info<vff::exceptions::integer>(e) )
          {
            std::ostringstream sstr;
            sstr << "Wrong site index: " << *mi << ".\n" << boost::diagnostic_information(e);
            PyErr_SetString(PyExc_ValueError, sstr.str().c_str());
          }
          else if(std::string const *mi = boost::get_error_info<vff::exceptions::string>(e))
            PyErr_SetString(PyExc_ValueError, mi->c_str());
          else PyErr_SetString(PyExc_ValueError, "Wrong site index in structure.");
          bp::throw_error_already_set();
        }
      }

    template<class T> 
      bp::tuple __call__( T &_self, LADA_MPI_CODE(boost::mpi::communicator const &_comm LADA_COMMA) bool relax )
      { 
        types::t_real const energy = _self.evaluate(LADA_MPI_CODE(_comm LADA_COMMA) relax) / 16.0217733;
        
        boost::shared_ptr< Crystal::TStructure<std::string> >
          result( new Crystal::TStructure<std::string>(_self.get_structure()) );
        result->energy = energy;
        if(relax) return bp::make_tuple(result, math::rMatrix3d(_self.get_stress()));
        else return bp::make_tuple(result, bp::object());
      }

    template<class T> 
      bp::tuple __call2__( T &_self, typename Crystal::TStructure<std::string> const &_str,
                           LADA_MPI_CODE(boost::mpi::communicator const &_comm LADA_COMMA)
                           bool _doinit, bool relax )
      { 
        init(_self, _str, _doinit);
        return __call__<T>(_self, LADA_MPI_CODE(_comm LADA_COMMA) relax);
      }

    template<class T> std::string angle_index(T &_self, bp::object const &_object)
    {
      if(bp::len(_object) != 3)
      {
        PyErr_SetString(PyExc_IndexError, "Angle parameters are accessed with 3-tuples.");
        bp::throw_error_already_set();
        return "";
      }
      std::string A, B, C;
      try
      {
        A = bp::extract<std::string>(_object[0]);
        B = bp::extract<std::string>(_object[1]);
        C = bp::extract<std::string>(_object[2]);
      }
      catch(...)
      {
        PyErr_SetString(PyExc_ValueError, "Could not convert index to 3-tuple of strings.");
        bp::throw_error_already_set();
        return "";
      }
      try { return vff::angle_type(B, A, C); }
      catch(vff::exceptions::angle_input &e)
      {
        PyErr_SetString(PyExc_ValueError, "Empty string is not a valid index.");
        bp::throw_error_already_set();
      }
      return "";
    }

    template<class T> std::string bond_index(T &_self, bp::object const &_object)
    {
      if(bp::len(_object) != 2)
      {
        PyErr_SetString(PyExc_IndexError, "Angle parameters are accessed with 3-tuples.");
        bp::throw_error_already_set();
        return "";
      }
      std::string A, B, index;
      try
      {
        A = bp::extract<std::string>(_object[0]);
        B = bp::extract<std::string>(_object[1]);
      }
      catch(...)
      {
        PyErr_SetString(PyExc_ValueError, "Could not convert index to 3-tuple of strings.");
        bp::throw_error_already_set();
        return "";
      }
      try { return vff::bond_type(A, B); }
      catch(vff::exceptions::bond_input &e)
      {
        PyErr_SetString(PyExc_ValueError, "Empty string is not a valid index.");
        bp::throw_error_already_set();
      }
      catch(std::exception &e)
      {
        PyErr_SetString(PyExc_RuntimeError, "Unknown error.");
        bp::throw_error_already_set();
      }
      return "";
    }

    template<class T> vff::BondData& get_bond(T &_self, bp::object const &_object)
    {
      std::string const index = bond_index(_self, _object);
      return _self.VffBase().bonds_params[index];
    }
    template<class T> void set_bond(T &_self, bp::object const &_index, vff::BondData const &_bond)
    {
      std::string const index = bond_index(_self, _index);
      try { _self.VffBase().bonds_params[index] = _bond; }
      catch(std::exception &e)
      {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        bp::throw_error_already_set();
      }
    }

    template<class T> vff::AngleData& get_angle(T &_self, bp::object const &_object)
    {
      std::string const index = angle_index(_self, _object);
      return _self.VffBase().angles_params[index];
    }
    template<class T> void set_angle(T &_self, bp::object const &_index, vff::AngleData const &_angle)
    {
      std::string const index = bond_index(_self, _index);
      try { _self.VffBase().angles_params[index] = _angle; }
      catch(std::exception &e)
      {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        bp::throw_error_already_set();
      }
    }
    template<class T> void check_input(T const &_self)
    {
      try { _self.VffBase().check_input(); }
      catch(vff::exceptions::missing_bond &e)
      {
        using vff::exceptions::bond;
        if(bond::value_type *ptr_bond = boost::get_error_info<bond>(e))
        {
          std::string const error = "Missing bond parameters for " + ptr_bond->get<0>() + "-" 
                                    + ptr_bond->get<1>() + " in vff input.";
          PyErr_SetString(PyExc_ValueError, error.c_str());
        }
        else PyErr_SetString(PyExc_ValueError, "Missing bond parameters in vff input.");
        bp::throw_error_already_set();
      }
      catch(vff::exceptions::missing_angle &e)
      {
        using vff::exceptions::angle;
        if(angle::value_type *ptr_angle = boost::get_error_info<angle>(e))
        {
          std::string const error = "Missing angle parameters for " + ptr_angle->get<1>() + "-" 
                                    + ptr_angle->get<0>() + "-" + ptr_angle->get<2>()
                                    + " in vff input.";
          PyErr_SetString(PyExc_ValueError, error.c_str());
        }
        else PyErr_SetString(PyExc_ValueError, "Missing angle parameters in vff input.");
        bp::throw_error_already_set();
      }
      catch(vff::exceptions::faulty_structure &e)
      {
        using vff::exceptions::atom;
        std::ostringstream sstr;
        if(atom::value_type *ptr_atom = boost::get_error_info<atom>(e))
          sstr << "Following atom is not four-fold coordinated.\n"
               << "   position:  (" << ptr_atom->pos[0] << ", "
               <<                      ptr_atom->pos[1] << ", "
               <<                      ptr_atom->pos[2] << "), type: "
               << ptr_atom->type << ".\n";
        else sstr << "Found atom which is not four-fold coordinated in structure.\n";
        sstr << "You may want to adjust bond_cutoff in vff, or correct the structure.\n";
        PyErr_SetString(PyExc_ValueError, sstr.str().c_str());
        bp::throw_error_already_set();
      }
    }

    template<class T> bp::class_<T> expose_functional(std::string const &_name, std::string const &_doc)
    {
      return bp::class_<T>(_name.c_str(), _doc.c_str())
        .def(bp::init<T const &>())
        .def_readwrite("_minimizer", &T::minimizer)
        .def( "_get_bond", &get_bond<T>, bp::return_internal_reference<>())
        .def( "_set_bond", &set_bond<T>)
        .def( "_get_angle", &get_angle<T>, bp::return_internal_reference<>())
        .def( "_set_angle", &set_angle<T>)
        .def
        ( 
          "__call__",  
          &__call__<T>,
          ( LADA_MPI_CODE(bp::arg("comm") LADA_COMMA) bp::arg("relax")=true )
        )
        .def
        ( 
          "__call__",  
          &__call2__<T>,
          ( bp::arg("structure"), LADA_MPI_CODE(bp::arg("comm") LADA_COMMA)
            bp::arg("doinit") = true, bp::arg("relax")=true ),
          "Minimizes structure.\n\n"
          ":Parameters:\n"
          "  structure : `crystal.Structure`\n"
          "    structure to evaluate.\n"
          LADA_MPI_CODE("  comm : boost.mpi.communicator\n    MPI communicator.\n")
          "  doinit : bool\n"
          "   If true, reconstructs first-neighbor tree. Default: true.\n\n"
          ":return: A 2-tuple (structure, stress) where the first is the relaxed "
          "structure, and the second a matrix with the stress. The energy is in "
          "``structure.energy``.\n"
        )
       .def("check_input", &check_input<T>)
       .def( "init",  &init<T>, (bp::arg("structure"), bp::arg("dotree") = true,
                                 bp::arg("verbose")=false), 
             "Initializes the functional for the current structure." ) 
       .def( "print_escan_input",  &print_escan_input<T>,
             (bp::arg("file"), bp::arg("structure") = bp::object()), 
             "Outputs the current structure in a format suitable for pescan." );
    }


    typedef vff::WithMinimizer< vff::Functional > t_Vff;
    typedef vff::WithMinimizer< vff::Layered > t_LayeredVff;
    math::rVector3d get_direction(t_LayeredVff const &_self)
      { return _self.VffBase().get_direction(); } 
    void set_direction(t_LayeredVff &_self, math::rVector3d const &_dir)
      { _self.VffBase().set_direction(_dir); } 

    
    void expose_vff()
    {
      expose_functional<t_Vff>
      ( 
        "Vff", 
        "A Valence Force Field Functional.\n\n"
        "Usage:\n"
        ">>> vff = lada.vff.Vff(\"input.xml\", boost.mpi.world)\n"
        ">>> structure = crystal.Structure(\"input.xml\")\n"
        ">>> relaxed_structure = vff(structure)\n\n"
      );
      bp::register_ptr_to_python< boost::shared_ptr<t_Vff> >();
    }

    void expose_layeredvff()
    {
      expose_functional<t_LayeredVff>
      ( 
        "LayeredVff", 
        "A Valence Force Field Functional with epitaxial constraints.\n\n"
        "Usage:\n"
        ">>> vff = lada.vff.LayeredVff(\"input.xml\", boost.mpi.world)\n"
        ">>> vff.direction = numpy.array([0,0,1], dtype=\"float64\")"
        ">>> structure = crystal.Structure(\"input.xml\")\n"
        ">>> relaxed_structure = vff(structure)\n\n"
      ).add_property
       (
         "direction",
         bp::make_function
         (
           &get_direction, 
           bp::return_value_policy<bp::return_by_value>()
         ),
         &set_direction, "Growth/Epitaxial direction.\n\n3x1 float64 numpy array.\n" 
       ); 
      bp::register_ptr_to_python< boost::shared_ptr<t_LayeredVff> >();
    }

  }
} // namespace LaDa
