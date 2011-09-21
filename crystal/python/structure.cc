#include "LaDaConfig.h"

#include <boost/python/class.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/make_function.hpp>
#include <boost/python/raw_function.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/with_custodian_and_ward.hpp>

#include <python/numpy_types.h>
#include <python/misc.hpp>
#include <python/exceptions.h>

#include "../structure.h"

#include "structure.hpp"
#include "add_atom.hpp"
#include "structure_contains.hpp"


namespace LaDa
{
  namespace python
  {
    namespace bp = boost::python;

    template<class T>
      bp::object get_cell(crystal::TemplateStructure<T> &_str ) 
        { return math::numpy::array_from_ptr<math::rMatrix3d::Scalar>(_str.cell().data(), 3u, 3u); }
    template< class T>
      void set_cell(crystal::TemplateStructure<T> &_str, bp::object _cell)
      {
        bp::object i_row( bp::handle<>(PyObject_GetIter(_cell.ptr())) );
        if(i_row.ptr() == Py_None)
        {
          PyErr_Clear();
          PyException<error::ValueError>::throw_error("Input is not a sequence.");
        }
        size_t i(0);
        for(; i < 3; ++i)
          if(PyObject *row = PyIter_Next(i_row.ptr()))
          {
            bp::object i_element( bp::handle<>(PyObject_GetIter(row)) );
            if(i_element.ptr() == Py_None)
            {
              PyErr_Clear();
              PyException<error::ValueError>::throw_error("Input is not a 3x3 sequence.");
            }
            size_t j(0);
            for(; j < 3; ++j)
              if(PyObject *element = PyIter_Next(i_element.ptr()))
              {
                if(PyFloat_Check(element))
                  _str.cell(i,j) = PyFloat_AS_DOUBLE(element);
                else if(PyInt_Check(element))
                  _str.cell(i,j) = (math::rMatrix3d::Scalar) PyInt_AS_LONG(element);
                else
                  PyException<error::ValueError>::throw_error("Input is not a 3x3 sequence of floats.");
                Py_DECREF(element);
              }
              else break;
            if(PyErr_Occurred()) return;
            if(j != 3) PyException<error::ValueError>::throw_error("Input is not a 3 by 3 sequence.");
            if(PyObject *element = PyIter_Next(i_element.ptr()))
            {
              Py_DECREF(element);
              PyException<error::ValueError>::throw_error
                ( "Input seems to be a 3xn sequence, with n > 3. Expected n == 3." );
            }
            Py_DECREF(row);
          }
          else break;
        if(PyErr_Occurred()) return;
        if(i != 3) PyException<error::ValueError>::throw_error("Input is not a 3 by 3 sequence.");
        if(PyObject *row = PyIter_Next(i_row.ptr()))
        {
          Py_DECREF(row);
          PyException<error::ValueError>::throw_error
            ( "Input seems to be a nx3 sequence, with n > 3. Expected n == 3." );
        }
      }
    template<class T> 
      types::t_real get_scale(crystal::TemplateStructure<T> &_in) { return _in.scale(); }
    template<class T> 
      void set_scale(crystal::TemplateStructure<T> &_in, types::t_real &_e) { _in.scale() = _e; }
    template<class T> 
      types::t_real get_energy(crystal::TemplateStructure<T> &_in) { return _in.energy(); }
    template<class T> 
      void set_energy(crystal::TemplateStructure<T> &_in, types::t_real &_e) { _in.energy() = _e; }
    template<class T> 
      types::t_real get_weight(crystal::TemplateStructure<T> &_in) { return _in.weight(); }
    template<class T> 
      void set_weight(crystal::TemplateStructure<T> &_in, bp::object const _e)
        { _in.weight() = bp::extract<types::t_real>(_e); }
    template<class T> 
      std::string get_name(crystal::TemplateStructure<T> &_in) { return _in.name(); }
    template<class T> 
      void set_name(crystal::TemplateStructure<T> &_in, std::string const &_e) { _in.name() = _e; }
    template<class T> 
      types::t_unsigned get_freeze(crystal::TemplateStructure<T> &_in) { return _in.freeze(); }
    template<class T> 
      void set_freeze(crystal::TemplateStructure<T> &_in, types::t_unsigned const &_e) { _in.freeze() = _e; }
    template<class T>
      std::string scalar_kind(crystal::TemplateStructure<T> const &) { return "scalar"; }
    template<class T>
      std::string list_kind(crystal::TemplateStructure<T> const &) { return "list"; }
    template<class T>
      std::string set_kind(crystal::TemplateStructure<T> const &) { return "set"; }

    template<class T>
      int __len__(crystal::TemplateStructure<T> const &_str) { return _str.size(); }
    template<class T>
      typename crystal::TemplateStructure<T>::reference 
        __getitem__(crystal::TemplateStructure<T> &_str, int _i)
        {
          if(_i < 0) _i += _str.size();
          if(_i < 0 or _i >= _str.size())
            PyException<error::IndexError>::throw_error("Index out of range.");
          return _str[_i];
        }
    //! Constructs and calls AddAtom functor.
    template<class T>
      bp::object add_atom(bp::tuple const &_args, bp::dict const &_kwargs)
      {
        crystal::TemplateStructure<T> structure = bp::extract< crystal::TemplateStructure<T> >(_args[0]);
        AddAtom<T> addatom(structure);
        bp::object result(addatom);
        addatom(bp::tuple(_args.slice(1, bp::slice_nil())), _kwargs);
        return bp::object(result);
      }

    template<class T>
      bp::class_< crystal::TemplateStructure<T> > 
        expose( std::string const &_name, std::string const &_desc ) 
        {
          return bp::class_<crystal::TemplateStructure<T> >( _name.c_str(), _desc.c_str() )
            .add_property
            (
              "cell",
              bp::make_function(&get_cell<T>, bp::with_custodian_and_ward_postcall<1, 0>()),
              &set_cell<T>, 
              "Cell vectors in cartesian coordinates.\n\n"
              "Units are ``self.scale``. The latter should be in angstrom, "
              "though use quantities package is not possible here."
            )
            .add_property( "energy",  &get_energy<T>, &set_energy<T>, "Energy of the structure." )
            .add_property( "scale",  &get_scale<T>, &set_scale<T>, "Scale of the cartesian units in Angstrom." )
            .add_property( "weight",  &get_weight<T>, &set_weight<T>, "Weight of the structure in fitting algorithms." )
            .add_property( "name",  &get_name<T>, &set_name<T>, "Name of the structure." )
            .add_property( "freeze", &get_freeze<T>, &set_freeze<T>,
                             "Tags to freeze coordinates when relaxing structure.\n\n" 
                             "See `FreezeCell` for possible values." 
                          )
            .def( "__len__", &__len__<T>, "Number of atoms in structure." )
            .def( "__getitem__", &__getitem__<T>, 
                  bp::return_internal_reference<>() )
            .def( "add_atom", bp::raw_function(&add_atom<T>, 1),
                  "Adds atom to structure.\n\n"
                  "See `lada.crystal.Atom` for parameter details."
                  "The ``kind`` parameter is not supported however (or meaningful). "
                  "The calls to this function can be chained as follows:\n\n"
                  ">>> structure.add_atom(0,0,0,\"Au\")\n"
                  ">>>                   (0.25,0.25,0.25,\"Pd\")\n" )
            .def( "__contains__", &details::contains<T>, 
                  "Returns True if a structure contains the given object.\n\n"
                  "The exact behavior depends on the type of the object and on the structure kind.\n"
                  " - list or set of species: ``structure.kind`` must 'list' or 'set'. "
                  "Returns True if an atomic site contains exactly that set of species.\n"
                  " - str: ``structure.kind`` must be 'scalar'. "
                  "Returns True if an atomic site is occupied by that specie.\n" 
                  " - sequence of three numbers (numpy array, list, ...): "
                  "Returns True if there exists a corresponding atomic site or periodic image. "
                  "In units of ``structure.scale``.\n"
                  " - object with 'pos' and 'type' attributes: combines the above.")
            .def_pickle( python::pickle< crystal::TemplateStructure<T> >() );
        }

    void expose_structure()
    {
      import_array(); // needed for NumPy 
      bp::enum_<crystal::frozenstr::type>( "FreezeCell", "Tags to freeze cell coordinates." )
        .value( "none", crystal::frozenstr::NONE )
        .value(   "xx", crystal::frozenstr::XX )
        .value(   "xy", crystal::frozenstr::XY )
        .value(   "xz", crystal::frozenstr::XZ )
        .value(   "yx", crystal::frozenstr::YX )
        .value(   "yy", crystal::frozenstr::YY )
        .value(   "yz", crystal::frozenstr::YZ )
        .value(   "zx", crystal::frozenstr::ZX )
        .value(   "zy", crystal::frozenstr::ZY )
        .value(   "zz", crystal::frozenstr::ZZ )
        .value(  "all", crystal::frozenstr::ALL )
        .value(  "a0", crystal::frozenstr::A0 )
        .value(  "a1", crystal::frozenstr::A1 )
        .value(  "a2", crystal::frozenstr::A2 )
        .export_values();

      expose_add_atom<std::string>("_AddAtomStr");
      expose<std::string>
      (
        "StructureStr", 
        "Defines a structure for wich atomic occupations are strings."
      ).add_property("kind", &scalar_kind<std::string>, "Occupations are strings.");
      expose_add_atom< std::vector<std::string> >("_AddAtomVec");
      expose< std::vector<std::string> >
      (
        "StructureVec", 
        "Defines a structure for wich atomic occupations are lists of strings."
      ).add_property("kind", &list_kind< std::vector<std::string> >, "Occupations is lists of strings.");
      expose_add_atom< std::set<std::string> >("_AddAtomSet");
      expose< std::set<std::string> >
      (
        "StructureSet", 
        "Defines a structure for wich atomic occupations are lists of strings."
      ).add_property("kind", &set_kind< std::set<std::string> >, "Occupations is sets of strings.");
    }

  }
} // namespace LaDa
