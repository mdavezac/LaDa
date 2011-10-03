#ifndef LADA_CRYSTAL_PYTHON_STRUCTURE_GETITEM_HPP
#define LADA_CRYSTAL_PYTHON_STRUCTURE_GETITEM_HPP

#include "LaDaConfig.h"

#include <math/python/python.hpp>

namespace LaDa
{
  namespace python
  {
    namespace details
    {
      template<class T> class ElementProxy;
    }
  }
}

namespace boost
{
  // Declared here in case compiler does not implement ADL.
  template <class T> inline typename LaDa::crystal::TemplateStructure<T>::value_type* 
    get_pointer(LaDa::python::details::ElementProxy<T> const& _p)
      { return _p.get(); }
}

namespace LaDa
{
  namespace python
  {
    namespace details
    {

      //! Proxy to atoms in a structure.
      template<class T> class ElementProxy
      {
        public:
          //! Type of the element this proxy points to.
          typedef typename crystal::TemplateStructure<T>::value_type element_type;
          //! Proxy constructor.
          ElementProxy   (crystal::TemplateStructure<T> const &_structure, size_t _i)
                       : structure_(_structure), i_(_i) {}
          //! Returns reference to element.
          element_type& operator*() const
          {
            if(i_ >= structure_.size())
              PyException<error::IndexError>::throw_error("Index-error when getting atom from structure.");
            return (*structure_.get())[i_];
          }
          //! Returns pointer to element.
          element_type* get() const
          {
            if(i_ >= structure_.size())
              PyException<error::IndexError>::throw_error("Index-error when getting atom from structure.");
            return &(*structure_.get())[i_];
          }
        private:
          //! Reference to the structure, keeping it alive.
          crystal::TemplateStructure<T> structure_;
          //! Index into the structure.
          size_t i_;
      };

  // Declared here in case compiler does not implement ADL.
  template <class T> inline typename LaDa::crystal::TemplateStructure<T>::value_type* 
    get_pointer(LaDa::python::details::ElementProxy<T> const& _p)
      { return _p.get(); }

      //! Suite defining all vector-like functions of a structure's python interface.
      template<class T> class StructureIndexing : public bp::def_visitor< StructureIndexing<T> >
      {
        public:
          //! Type of the structure being looked at.
          typedef crystal::TemplateStructure<T> t_Structure;
          //! Visitor where additional methods are added to the class.
          template<class CLASS> void visit(CLASS &_cl) const
          {
            bp::register_ptr_to_python< ElementProxy<T> >();
            _cl.def("__len__", &getitem)
               .def("__getitem__", &getitem)
               .def("__delitem__", &delitem)
               .def("__contains__", &contains);
          };

          //! Returns size of container.
          static int len(t_Structure const &_structure) { return _structure.size(); }
          //! Returns item or slice.
          static bp::object getitem(t_Structure & _structure, PyObject *_i) 
          {
            bp::object const i(bp::borrowed(_i));
            if(PyInt_Check(_i))
            {
              int i = PyInt_AS_LONG(_i);
              if(i < 0) i += _structure.size();
              if(i < 0 or i >= _structure.size())
                PyException<error::IndexError>::throw_error("Index out of range in Structure's __getitem__.");
              return bp::object(ElementProxy<T>(_structure, i));
            } 
            else if(is_position(i))
            {
              typename t_Structure::iterator const i_found = getitem_pos_(_structure, i); 
              if(i_found == _structure.end()) 
                PyException<error::IndexError>::throw_error("Could not find corresponding atomic position.");
              return bp::object(ElementProxy<T>(_structure, i_found - _structure.begin()));
            }
            PyException<error::TypeError>::throw_error("Unknown argument to Structure's __getitem__.");
          }
          //! Deletes item or slice.
          static void delitem(t_Structure & _structure, PyObject *_i) 
          {
            bp::object i(bp::borrowed<>(_i));
            if(PyInt_Check(_i))
            {
              int i = PyInt_AS_LONG(_i);
              if(i < 0) i += _structure.size();
              if(i < 0 or i >= _structure.size())
                PyException<error::IndexError>::throw_error("Index out of range in Structure's __getitem__.");
              _structure.erase(_structure.begin()+i);
              return;
            } 
            else if(is_position(i))
            {
              typename t_Structure::iterator const i_found = getitem_pos_(_structure, i);
              if(i_found == _structure.end()) 
                PyException<error::IndexError>::throw_error("Could not find corresponding atomic position.");
              _structure.erase(i_found);
              return;
            }
            else if(is_specie_index<T>(i))
            {
              T filter;
              extract_species(i, filter);
              crystal::CompareOccupations<T> const cmp(filter);
              typename crystal::TemplateStructure<T>::reverse_iterator i_atom = _structure.rbegin();
              typename crystal::TemplateStructure<T>::reverse_iterator const i_atom_end = _structure.rend();
              bool found_one(false);
              for(size_t n(_structure.size()-1); i_atom != i_atom_end; ++i_atom, --n)
                if(cmp(i_atom->type))
                {
                  found_one = true;
                  _structure.erase(_structure.begin() + n);
                }
              if(not found_one)
                PyException<error::IndexError>::throw_error("Could not find corresponding specie(s).");
              return;
            }
            PyException<error::TypeError>::throw_error("Unknown argument to Structure's __delitem__.");
          }
          //! True if item is contained in structure.
          static bool contains(t_Structure & _structure, PyObject *_i) 
          {
            bp::object const i(bp::borrowed(_i));
            if(is_position(bp::object(bp::borrowed<>(_i))))
              return _structure.end() != getitem_pos_(_structure, i);
            else if(is_specie_index<T>(i))
            {
              T filter;
              extract_species(i, filter);
              crystal::CompareOccupations<T> const cmp(filter);
              typename crystal::TemplateStructure<T>::const_iterator i_atom = _structure.begin();
              typename crystal::TemplateStructure<T>::const_iterator const i_atom_end = _structure.end();
              for(; i_atom != i_atom_end; ++i_atom)
                if(cmp(i_atom->type)) return true;
              return false;
            }
            PyException<error::TypeError>::throw_error("Unknown argument to Structure's __contains__.");
          }
        private:
          //! Returns atom from knowledge of position.
          static typename t_Structure::iterator
            getitem_pos_(t_Structure &_structure, bp::object const &_index)
            {
              bp::extract<math::rVector3d::Scalar const&> const x(_index[0]),
                                                                y(_index[1]),
                                                                z(_index[2]);
              math::rVector3d const filter(x(), y(), z());
              typename crystal::TemplateStructure<T>::iterator i_atom = _structure.begin();
              typename crystal::TemplateStructure<T>::iterator const i_atom_end = _structure.end();
              math::rMatrix3d const inv(_structure.cell().inverse());
              for(; i_atom != i_atom_end; ++i_atom)
                if( math::is_integer(inv * (i_atom->pos - filter)) ) return i_atom;
              PyException<error::IndexError>::throw_error("No atom found at given position.");
              return i_atom_end;
            }
      };

    }
  }
}
#endif
