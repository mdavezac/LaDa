#ifndef LADA_CRYSTAL_ADD_ATOM_HPP
#define LADA_CRYSTAL_ADD_ATOM_HPP

#include "LaDaConfig.h"

#include <boost/python/return_arg.hpp> 

#include "atom.hpp"

namespace LaDa
{
  namespace python
  {
    //! Functor to add an atom to a structure.
    template<class T>
      struct AddAtom
      {
        //! Reference to structure to which atoms are added.
        crystal::TemplateStructure<T> structure;
        //! \brief Constructor.
        //! \details Should only be called from c++. 
        AddAtom(crystal::TemplateStructure<T> &_str) : structure(_str) {}
        //! Copy Constructor.
        AddAtom(AddAtom<T> const &_c) : structure(_c.structure) {}
        //! \brief Actually adds atom to structure.
        //! \details Should follow `lada.crystal.Atom` input parameters.
        void operator()(bp::tuple const &_args, bp::dict const &_kwargs);
      };

    template<class T>
      void AddAtom<T>::operator()(bp::tuple const &_args, bp::dict const &_kwargs) 
      {
        typedef typename crystal::TemplateStructure<T>::value_type t_Atom;
        if(bp::len(_args) == 0 and bp::len(_kwargs) == 0) 
          structure.push_back(t_Atom());
        t_Atom atom;
        if(bp::len(_args) >= 3 and _kwargs.has_key("position"))
          PyException<error::TypeError>::throw_error
             ("Position provided both in tuple and keyword argument.");
        else if(bp::len(_args) >= 3 or _kwargs.has_key("position"))
        {
          bp::object pos;
          if(_kwargs.has_key("position")) pos = _kwargs.get("position");
          else pos = _args.slice(bp::slice_nil(), 3);
          if(bp::len(pos) != 3) 
            PyException<error::TypeError>::throw_error
              ("Position is not a list of three scalars.");
          atom.pos[0] = boost::python::extract<math::rVector3d::Scalar>(pos[0]);
          atom.pos[1] = boost::python::extract<math::rVector3d::Scalar>(pos[1]);
          atom.pos[2] = boost::python::extract<math::rVector3d::Scalar>(pos[2]);
        }
        if(bp::len(_args) > 3 and _kwargs.has_key("type"))
          PyException<error::TypeError>::throw_error
             ("Types provided both in tuple and keyword argument.");
        else if(bp::len(_args) == 4 or _kwargs.has_key("type"))
        {
          bp::object type;
          if(_kwargs.has_key("type")) type = _kwargs.get("type");
          else type = _args[3];
          details::extract_type0(type, atom.type);
        }
        else if(bp::len(_args) > 4)
          details::extract_type1(_args.slice(3, bp::slice_nil()), atom.type);
        if(_kwargs.has_key("freeze"))
          atom.freeze = bp::extract<int>(_kwargs.get("freeze"));
        if(_kwargs.has_key("site"))
          atom.site = bp::extract<int>(_kwargs.get("site"));
        structure.push_back(atom);
      }

    template<class T>
      bp::object add_atom_(bp::tuple const &_args, bp::dict const &_kwargs)
      {
        bp::object self(_args[0]);
        AddAtom<T> addatom = bp::extract< AddAtom<T> >(self);
        addatom(bp::tuple(_args.slice(1, bp::slice_nil())), _kwargs);
        return self;
      }

    template<class T> void expose_add_atom(std::string const _name)
    {
      bp::class_< AddAtom<T> >(_name.c_str(), "Class to add atoms to a structure.", bp::no_init)
        .def( "__call__", bp::raw_function(&add_atom_<T>),
              "Adds atom to structure.\n\n"
              "See `lada.crystal.Atom` for parameter details."
              "The ``kind`` parameter is not supported however (or meaningful)." );
    }
  }
}

#endif 
