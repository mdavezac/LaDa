#ifndef LADA_CRYSTAL_PYTHON_STRUCTURE_CONTAINS_HPP
#define LADA_CRYSTAL_PYTHON_STRUCTURE_CONTAINS_HPP

#include "LaDaConfig.h"

#include <iterator> 

#include <boost/utility/enable_if.hpp>
#include <boost/python/stl_iterator.hpp>

#include <math/misc.h>
#include "../is_container.h"
#include "../compare_sites.h"

namespace LaDa
{
  namespace python
  {
    namespace details
    {
      // extract sequence to specie.
      template<class T> typename boost::enable_if< crystal::details::is_scalar<T>, bool > :: type
        extract_species(bp::object const &_filter, T& _in)
          { _in = bp::extract<T>(_filter); }
      // extract sequence to specie.
      template<class T> typename boost::enable_if< crystal::details::is_container<T>, bool > :: type
        extract_species(bp::object const &_filter, T& _in)
        {
          std::back_insert_iterator<T> i_input(_in);
          bp::stl_input_iterator<typename T::value_type> i_filter(_filter);
          bp::stl_input_iterator<typename T::value_type> const i_filter_end;
          std::copy(i_filter, i_filter_end, i_input);
        }
      // extract sequence to specie.
      template<class T> typename boost::enable_if< crystal::details::is_set<T>, bool > :: type
        extract_species(bp::object const &_filter, T& _in)
        {
          bp::stl_input_iterator<typename T::value_type> i_filter(_filter);
          bp::stl_input_iterator<typename T::value_type> const i_filter_end;
          for(; i_filter != i_filter_end; ++i_filter) _in.insert(*i_filter);
        }

      void extract_position(bp::object const &_in, math::rVector3d &_out)
      {
        bp::stl_input_iterator<math::rVector3d::Scalar> i_filter(_in);
        bp::stl_input_iterator<math::rVector3d::Scalar> const i_filter_end;
        for(size_t i(0); i < 3; ++i, ++i_filter)
        { 
          if(i_filter == i_filter_end)
            PyException<error::TypeError>::throw_error("Not a sequence of three numbers.");
          _out[i] = *i_filter;
        }
      }

      //! Checks if a specie is in a list structure.
      template<class T> 
        bool __contains_specie__(crystal::TemplateStructure<T> const &_in, bp::object const &_filter)
        {
          T filter;
          extract_species(_filter, filter);
          crystal::CompareOccupations<T> const cmp(filter);
          typename crystal::TemplateStructure<T>::const_iterator i_atom = _in.begin();
          typename crystal::TemplateStructure<T>::const_iterator const i_atom_end = _in.end();
          for(; i_atom != i_atom_end; ++i_atom)
            if(cmp(i_atom->type)) return true;
          return false;
        }
      //! Checks if a position is in a list structure.
      template<class T> 
        bool __contains_pos__(crystal::TemplateStructure<T> const &_in, bp::object const &_filter)
        {
          math::rVector3d filter;
          extract_position(_filter, filter);

          typename crystal::TemplateStructure<T>::const_iterator i_atom = _in.begin();
          typename crystal::TemplateStructure<T>::const_iterator const i_atom_end = _in.end();
          math::rMatrix3d const inv(_in.cell().inverse());
          for(; i_atom != i_atom_end; ++i_atom)
            if( math::is_integer(inv * (i_atom->pos - filter)) ) return true;
          return false;
        }

      //! Cannot check that an empty list of species is in a scalar structure.
      template<class T> typename boost::enable_if< crystal::details::is_scalar<T>, bool > :: type
        __contains_empty__(crystal::TemplateStructure<T> const &)
        {
          PyException<error::ValueError>::throw_error
            ("Cannot use ``a in structure`` if a is a list and structure scalar.");
          return false;
        }
      //! Checks if a list of species is in a list structure.
      template<class T> typename boost::enable_if< crystal::details::is_iterable<T>, bool > :: type
        __contains_empty__(crystal::TemplateStructure<T> const &_in)
        {
          typename crystal::TemplateStructure<T>::const_iterator i_atom = _in.begin();
          typename crystal::TemplateStructure<T>::const_iterator const i_atom_end = _in.end();
          for(; i_atom != i_atom_end; ++i_atom)
            if(i_atom->type.size() == 0) return true;
          return false;
        }
      //! Checks if an atom is contained in a structure.
      template<class T> 
        bool __contains_atom__(crystal::TemplateStructure<T> const &_in, bp::object const &_filter)
        {
          T specie;
          bp::object const sp_attr = _filter.attr("type");
          bp::object const pos_attr = _filter.attr("pos");
          extract_species(sp_attr, specie);
          math::rVector3d pos;
          extract_position(pos_attr, pos);

          typename crystal::TemplateStructure<T>::const_iterator i_atom = _in.begin();
          typename crystal::TemplateStructure<T>::const_iterator const i_atom_end = _in.end();
          math::rMatrix3d const inv(_in.cell().inverse());
          crystal::CompareOccupations<T> const cmp(specie);
          for(; i_atom != i_atom_end; ++i_atom)
            if( math::is_integer(inv * (i_atom->pos - pos)) 
                and cmp(i_atom->type) ) return true;
          return false;
        }




      //! Forks to different functions depending on input.
      template<class T>
        bool contains(crystal::TemplateStructure<T> const &_in, bp::object const &_filter)
        {
          if(     PyObject_HasAttrString(_filter.ptr(), "pos")
              and PyObject_HasAttrString(_filter.ptr(), "type") )
            return __contains_atom__<T>(_in, _filter);
          if(PySequence_Check(_filter.ptr()))
          {
            int const N(bp::len(_filter));
            if(N == 3) 
            {
              bp::object const item = _filter[0];
              if(PyInt_Check(item.ptr()) or PyFloat_Check(item.ptr()))
                return __contains_pos__<T>(_in, _filter);
            }
            else if(N == 0) return __contains_empty__<T>(_in);
          }
          return __contains_specie__<T>(_in, _filter);
        };

    }
  }
}
#endif
