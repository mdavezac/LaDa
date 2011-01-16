#include "LaDaConfig.h"


#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/data_members.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#ifdef LADA_MPI
# include <boost/mpi/collectives.hpp>
#endif


#include "vff.hpp"
#include "../functional.h"
#include "../layered.h"
#include "../va.h"

namespace LaDa
{
  namespace bp = boost::python;
  namespace python
  {
    template<class T> 
      bp::tuple __call__( T &_self, typename Crystal::TStructure<std::string> const &_str,
                          LADA_MPI_CODE(boost::mpi::communicator const &_comm LADA_COMMA)
                          bool _doinit, bool relax )
      { 
        _self.VffBase().structure = _str;
        _self.init(_doinit, false);
        types::t_real const energy = _self.evaluate(relax LADA_MPI_CODE(LADA_COMMA _comm)) / 16.0217733;
        
        boost::shared_ptr< Crystal::TStructure<std::string> >
          result( new Crystal::TStructure<std::string>(_self.structure) );
        result->energy = energy;
        if(relax) return bp::make_tuple(result, math::rMatrix3d(_self.first.get_stress()));
        else return bp::make_tuple(result, bp::object());
      }

    template<class T> std::map<std::string, vff::BondData> const & get_bonddata(T const &_t)
      {return _t.VffBase().bonds_params;}
    template<class T> void set_bonddata(T &_t, std::map<std::string, vff::BondData> const & _val)
      {_t.VffBase().bonds_params = _val;}
    template<class T> std::map<std::string, vff::AngleData> const & get_angledata(T const &_t)
      {return _t.VffBase().angles_params;}
    template<class T> void set_angledata(T &_t, std::map<std::string, vff::AngleData> const & _val)
      {_t.VffBase().angles_params = _val;}

    template<class T> bp::class_<T> expose_functional(std::string const &_name, std::string const &_doc)
    {
      return bp::class_<T>(_name.c_str(), _doc.c_str())
        .def(bp::init<T const &>())
        .def_readwrite("minimizer", &T::minimizer)
        .add_property( "_bonds_params", bp::make_function(&get_bonddata<T>, bp::return_internal_reference<>()), 
                       &set_bonddata<T> )  
        .add_property( "_angles_params", bp::make_function(&get_angledata<T>, bp::return_internal_reference<>()), 
                       &set_angledata<T> )  
        .def
        ( 
          "__call__",  
          &__call__<T>,
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
       .def( "_init",  &T::init, (bp::arg("redo_tree") = true, bp::arg("verbose")=false), 
             "Initializes the functional for the current structure." ) 
       .def( "print_escan_input",  &print_escan_input<T>,
             (bp::arg("file"), bp::arg("structure")), 
             "Outputs the current structure in a format suitable for pescan." );
    }


    typedef vff::VABase< vff::Functional > t_Vff;
    typedef vff::VABase< vff::Layered > t_LayeredVff;
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
