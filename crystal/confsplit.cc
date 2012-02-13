#include "LaDaConfig.h"

#include <algorithm>

#include <boost/bind.hpp>
#include <boost/ref.hpp>

#include <math/fuzzy.h>
#include <python/random_access_list_iterator.h>
#include <python/wrap_numpy.h>

#include "coordination_shells.h"
#include "confsplit.h"

namespace LaDa
{
  namespace crystal
  {
#   ifdef LADA_ADD_ITEM
#     error LADA_ADD_ITEM already defined
#   endif
#   define LADA_ADD_ITEM(BITSET, INDEX, OBJECT)  \
      if(PyList_SET_ITEM(BITSET.borrowed(), INDEX, OBJECT) != 0) return NULL;
#   ifdef LADA_GET_POS
#     error LADA_GET_POS already defined
#   endif
#   define LADA_GET_POS(NEIGH) ((AtomData*)PyTuple_GET_ITEM((PyTupleObject*)NEIGH, 0))->pos
    namespace 
    {
      //! creates new bitset.
      PyObject* new_bitset(Py_ssize_t _n)
      {
        PyObject *result = PyList_New(2);
        PyObject *list = PyList_New(_n);
        if(not list) goto error;
        if(PyList_SET_ITEM(result, 0, list) != 0) goto error;
        error:
          Py_DECREF(result);
          return NULL;
      }
      
      //! Look for largest x element.
      python::RAList_iterator max_xelement( python::RAList_iterator & _first,
                                            python::RAList_iterator const & _last,
                                            math::rVector3d const &_skip,
                                            math::rVector3d const &_x )
      {
        if (_first == _last) return _first;
        do
        {
          math::rVector3d const &vec(LADA_GET_POS(*_first));
          if( not math::is_null( (_skip - vec).squaredNorm()) ) continue;
        }
        while (++_first != _last);
        if(_first == _last) return _last;
        python::RAList_iterator _result = _first;
        while (++_first != _last)
        {
          math::rVector3d const &vec(LADA_GET_POS(*_first));
          if(math::is_null( (_skip - vec).squaredNorm())) continue;
          if(math::leq(LADA_GET_POS(*_result).dot(_x), vec.dot(_x))) _result = _first;
        }
        return _result;
      }
      
      //! Functor to compare coordinates using once given a basis.
      struct CmpFromCoord
      {
        math::rVector3d const &x;
        math::rVector3d const &y;
        math::rVector3d const &z;
        CmpFromCoord   (math::rVector3d const &_x, math::rVector3d const &_y, math::rVector3d const &_z) 
                     : x(_x), y(_y), z(_z) {}
        CmpFromCoord(CmpFromCoord const &_in) : x(_in.x), y(_in.y), z(_in.z) {}
        bool operator()(PyObject* const _a, PyObject* const _b) const
        {
          math::rVector3d const &a = ((AtomData*)PyTuple_GET_ITEM(_a, 0))->pos;
          math::rVector3d const &b = ((AtomData*)PyTuple_GET_ITEM(_b, 0))->pos;
          const types::t_real x1(a.dot(x)), x2(b.dot(x));
          if( math::neq(x1, x2) ) return math::gt(x1, x2);
          const types::t_real y1(a.dot(y)), y2(b.dot(y));
          if( math::neq(y1, y2)) return math::gt(y1, y2);
          return math::gt(a.dot(z), b.dot(z) );
        }
      };
    } // unnamed namespace.

    // Function to compare configurations.
    bool compare_configurations(PyObject* const _a, PyObject* const _b, types::t_real _tolerance)
    {
      if(PyList_GET_SIZE(_a) != PyList_GET_SIZE(_b)) return false;
      python::RAList_iterator i_first(_a, 0);
      python::RAList_iterator i_cmp(_b, 0);
      python::RAList_iterator i_end(_a);
      for(; i_first != i_end; ++i_first, ++i_cmp)
      {
        if(PyList_GET_ITEM(*i_first, 0) != PyList_GET_ITEM(*i_cmp, 0)) return false;
        types::t_real const a = PyFloat_AS_DOUBLE(PyList_GET_ITEM(*i_first, 2));
        types::t_real const b = PyFloat_AS_DOUBLE(PyList_GET_ITEM(*i_cmp, 2));
        if( math::neq(a, b, _tolerance) ) return false;
        Eigen::Map<math::rVector3d> const
          veca((math::rVector3d::Scalar*)PyArray_DATA(PyList_GET_ITEM(*i_first, 1))),
          vecb((math::rVector3d::Scalar*)PyArray_DATA(PyList_GET_ITEM(*i_cmp, 1)));
        if(math::neq(veca[0], vecb[0], _tolerance)) return false;
        if(math::neq(veca[1], vecb[1], _tolerance)) return false;
        if(math::neq(veca[2], vecb[2], _tolerance)) return false;
      }
      return true;
    };

    bool splitconfigs( Structure const &_structure,
                       Atom const &_origin,
                       Py_ssize_t _nmax,
                       python::Object &_configurations,
                       types::t_real const _tolerance )
    {
      const types::t_real weight( 1e0 / types::t_real(_structure.size()) );

      // if configuration is null, creates it.
      if(not _configurations)
      {
        _configurations.reset(PyList_New(0));
        if(not _configurations) return false;
      }

      // Creates new bitset and set origin as first of its references.
      python::Object bitset = new_bitset(_nmax);
      if(bitset)
      {
        PyObject* pydist = PyFloat_FromDouble(0);
        if(not pydist) return false;
        math::rVector3d const zero = math::rVector3d::Zero();
        PyObject* pytrans = python::wrap_to_numpy(zero);
        if(not pytrans) {Py_DECREF(pydist); return false; }
        PyObject* dummy = PyTuple_Pack(3, _origin.borrowed(), pytrans, pydist);
        Py_DECREF(pydist); Py_DECREF(pytrans);
        if(not dummy) return false;
        LADA_ADD_ITEM(bitset, 0, dummy);
      }
      else return false;

      const math::rVector3d origin(_origin->pos);

      python::Object epositions = coordination_shells(_structure, _nmax, origin, _tolerance);

      // loop over epositions defining x.
      python::RAList_iterator i_xpositions(epositions.borrowed(), 0);
      const types::t_real
        xweight( weight / types::t_real(PyList_GET_SIZE(*i_xpositions)) );
      python::RAList_iterator i_xpos(*i_xpositions, 0);
      python::RAList_iterator i_xpos_end(*i_xpositions);
      for(; i_xpos != i_xpos_end; ++i_xpos)
      {
        math::rVector3d const &xpos(LADA_GET_POS(*i_xpos));
        math::rVector3d const x(xpos - origin);

        // finds positions defining y.
        // Stores possible y positions.
        std::vector<PyObject*> ypossibles;
        python::RAList_iterator i_ypositions = i_xpositions;
        if( PyList_GET_SIZE(*i_xpositions) == 1 ) ++i_ypositions; 

        python::RAList_iterator i_ypos = python::RAList_iterator(*i_ypositions);
        python::RAList_iterator const i_ypos_end = python::RAList_iterator(*i_ypositions);
        python::RAList_iterator max_x_element = max_xelement( i_ypos, i_ypos_end, x, xpos);
        if(max_x_element == i_ypos_end) 
        {
          LADA_PYERROR(InternalError, "Should not be here.");
          return false;
        }
        const types::t_real max_x_scalar_pos(LADA_GET_POS(*max_x_element).dot(x) );
        for(i_ypos = python::RAList_iterator(*i_ypositions, 0); i_ypos != i_ypos_end; ++i_ypos)
        {
          math::rVector3d const &ypos = LADA_GET_POS(*i_ypos);
          if( math::neq(ypos.dot(x), max_x_scalar_pos) ) continue;
          if( math::is_null( (ypos - xpos).squaredNorm() ) ) continue;
          ypossibles.push_back(*i_ypos);
        }

        // divide current weight by number of possible y positions.
        types::t_real const bitsetweight = xweight / types::t_real(ypossibles.size());
        std::vector<PyObject*>::const_iterator i_ypossible = ypossibles.begin();
        std::vector<PyObject*>::const_iterator const i_ypossible_end = ypossibles.end();
        // loop over possible ys.
        //   Determine z from x and y.
        //   Basis is determined. Adds other atoms.
        for(; i_ypossible != i_ypossible_end; ++i_ypossible)
        {
          // at this point, we can define the complete coordinate system.
          const math::rVector3d y( ((AtomData*)PyTuple_GET_ITEM(*i_ypossible, 0))->pos - origin );
          const math::rVector3d z( x.cross(y) );

          // atoms are now included in the list according to the following rule:
          //  _ closest to the origin first.
          //  _ ties are broken according to largest x coordinate.
          //  _ next ties are broken according to largest y coordinate.
          //  _ final ties are broken according to largest z coordinate.

          // we iterate over coordination shells and add reference until nmax is reached.
          Py_ssize_t current_index = 1;
          python::RAList_iterator i_shell(epositions.borrowed(), 0);
          python::RAList_iterator const i_shell_end(epositions.borrowed());
          for(; i_shell != i_shell_end; ++i_shell)
          {
            if( current_index == _nmax ) break;

            const size_t edn( std::min(PyList_GET_SIZE(*i_shell), _nmax - current_index) );
            if(edn == 1) // case where the shell contains only one atom ref.
            {
              LADA_ADD_ITEM(bitset, current_index++, *i_ypos);
              continue;
            }

            // case where all atom in shell should be added.
            if( edn <= _nmax - current_index ) 
              std::sort
              ( 
                python::RAList_iterator(*i_shell, 0), python::RAList_iterator(*i_shell), 
                CmpFromCoord(x, y, z) 
              );
            // case where only a few atoms in the shell should be added.
            else std::partial_sort
                 ( 
                   python::RAList_iterator(*i_shell, 0),
                   python::RAList_iterator(*i_shell, _nmax),
                   python::RAList_iterator(*i_shell), 
                   CmpFromCoord(x, y, z) 
                 );
            python::RAList_iterator i_addme(*i_shell, 0);
            python::RAList_iterator const i_addme_end(*i_shell, edn);
            for(; i_addme != i_addme_end; ++i_addme, ++current_index)
            {
              Py_INCREF(*i_addme);
              LADA_ADD_ITEM(bitset, current_index, *i_addme);
            }
          } // end of loop over positions at equivalent distance.


          // finally adds configuration.
          python::RAList_iterator const i_conf_end(_configurations.borrowed());
          python::RAList_iterator i_found = std::find_if
            (
              python::RAList_iterator(_configurations.borrowed(), 0), i_conf_end,
              boost::bind(compare_configurations, _1, bitset.borrowed(), _tolerance)
            );
          if(i_found == i_conf_end) // add new bitset.
          {
            PyObject *dummy = PyFloat_FromDouble(bitsetweight); 
            if(not dummy) return false;
            // should not have been set yet.
            if(PyList_SET_ITEM(bitset.borrowed(), 1, dummy) != 0) return false;
            PyList_Append(_configurations.borrowed(), bitset.borrowed());
          }
          else // add to pre-existing bitset.
          {
            double const real = PyFloat_AS_DOUBLE(PyList_GET_ITEM(*i_found, 1));
            if(PyErr_Occurred() != NULL) return false;
            PyObject *dummy = PyFloat_FromDouble(bitsetweight); 
            if(not dummy) return false;
            if(PyList_SET_ITEM(*i_found, 1, dummy) != 0) return false;
          }
        } // end of loop over equivalent y coords.
      } // end of loop over equivalent  x coords.
      return true;
    }
#   undef LADA_GET_POS
#   undef LADA_ADD_ITEM
  } // namespace Crystal
} // namespace LaDa
