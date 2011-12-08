#ifndef LADA_CRYSTAL_ATOM_BASE_H
#define LADA_CRYSTAL_ATOM_BASE_H

#include "LaDaConfig.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <ostream>
#include <complex>
#include <iostream>

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#ifdef LADA_DO_PYTHON
# include <Python.h>
#endif

#include <opt/types.h>
#include <math/eigen.h>
#include <math/fuzzy.h>

#include <math/serialize.h>
#ifdef LADA_WITH_LNS
#  include "load_n_save/xpr/utilities.h"
#  include "load_n_save/action/vector.h"
#  include "load_n_save/action/set.h"
#endif
#include "is_container.h"
#include "atom_freeze.h"

namespace LaDa
{
  namespace crystal
  {
    //!\cond
    template<class T> class AtomData;
    //!\endcon


    //! Dumps atom to stream.
    template<class T> std::ostream& operator<<(std::ostream&, AtomData<T>const &);

    //! \brief Describes an atom.
    //! \details An atom consists of a position and a type. The position should
    //!          always be in cartesian units. The type can be anything, from a
    //!          string with the symbol of the atom, to an double wich codes
    //!          for the atomic type somehow, to a vector of strings which
    //!          describe the possible occupations of the atomic position. To
    //!          this end, the type is a template type \a T_TYPE. 
    //! \warning The default equality comparison operator compares positions only (not
    //!          occupation or site ).
    template<class T_TYPE>
      class AtomData: public AtomFreezeMixin
      {
        friend class boost::serialization::access;
#       ifdef LADA_WITH_LNS
          friend class load_n_save::access;
#       endif
        template<class T> friend std::ostream& operator<<(std::ostream&, AtomData<T>const &);
        public:
          //! The type of the occupation
          typedef T_TYPE t_Type;
          //! The atomic position in cartesian coordinate.
          math::rVector3d pos;
          //! The atomic occupation.
          t_Type  type;
          //! \brief Site index.
          //! \details Used to reference atomic sites in a supercell versus
          //!          atomic sites in a reference lattice.
          types::t_int site;
#         ifdef LADA_PYDICT
#           error LADA_PYDICT already defined.
#         endif
#         ifdef LADA_DO_PYTHON
            //! \brief Reference to python attribute dictionary.
            //! \details A reference is owned by this object. This way, the
            //!          dictionary of attributes, once created is never lost
            //!          over the course of the atom's life, whatever may
            //!          happen to python wrappers around the atom.
            PyObject *pydict;
            //! \brief Borrowed reference to a python wrapper around this object.
            //! \details To avoid impossible ownership cycles, this reference
            //!          is not owned by the atom object. Rather, whenever a
            //!          wrapper is created, the reference is checked for
            //!          existence first. If it is NULL, a new wrapper is
            //!          created and a this reference is made to point to it.
            //!          However, it is not incref'ed and thus not owned. Upon
            //!          destruction of the python wrapper, this reference is
            //!          set to NULL. By carefully using only the provided
            //!          PyAtom_FromAtom, it is possible to ensure that there
            //!          are never more than one wrapper around, and that
            //!          pyself always points to a valid memory.
            PyObject *pyself;
#           define LADA_PYDICT , pydict(NULL), pyself(NULL)
#         else 
#           define LADA_PYDICT
#         endif

          //! Constructor
          AtomData() : AtomFreezeMixin(frozen::NONE), pos(math::rVector3d(0,0,0)),
                       type(), site(-1) LADA_PYDICT {}
          //! Constructor
          template<class T_DERIVED>
            AtomData   ( Eigen::DenseBase<T_DERIVED> const &_pos, t_Type const &_in,
                         types::t_int _site = -1, types::t_unsigned _freeze = frozen::NONE )
                     : AtomFreezeMixin(_freeze), pos(_pos), type(_in), site(_site) LADA_PYDICT {}
          //! Copy Constructor
          AtomData   (const AtomData &_c)
                   : AtomFreezeMixin(_c), pos(_c.pos), type(_c.type), site(_c.site)
          {
#           ifdef LADA_DO_PYTHON
              pyself = NULL;
              if(_c.pydict == NULL and pydict != NULL) 
              {
                PyObject *dummy = pydict;
                pydict = NULL;
                Py_XDECREF(dummy);
              }
              else if(_c.pydict != NULL)
              {
                if(pydict != NULL) { Py_DECREF(pydict); pydict = NULL; }
                PyObject* copymod = PyImport_ImportModule("copy");
                if(copymod == NULL) return;
                PyObject *deepcopystr = PyString_FromString("deepcopy");
                if(not deepcopystr) { Py_DECREF(copymod); return; }
                pydict = PyObject_CallMethodObjArgs(copymod, deepcopystr, _c.pydict, NULL);
                Py_DECREF(copymod);
                Py_DECREF(deepcopystr);
              }
#           endif
          }
#         undef LADA_PYDICT
#         ifdef LADA_DO_PYTHON
          ~AtomData()
          { 
            if(pydict != NULL)
            {
              PyObject *dummy = pydict;
              pydict = NULL;
              Py_DECREF(dummy); 
            }
          }
#         endif
      
        private:
          //! Serializes an atom.
          template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
          {
             namespace bs = boost::serialization;
             _ar & bs::base_object< AtomData<T_TYPE> >(*this);
             _ar & pos; _ar & type; _ar & site;
          }
#         ifdef LADA_WITH_LNS
            //! To load and save to xml-like input.
            template<class T_ARCHIVE> bool lns_access(T_ARCHIVE &_ar, const unsigned int _version);
#         endif
      };

#   ifdef LADA_WITH_LNS
      //! To load and save to xml-like input.
      template<class T_TYPE> template<class T_ARCHIVE>
        bool AtomData<T_TYPE> :: lns_access(T_ARCHIVE &_ar, const unsigned int _version) 
        {
          namespace lns = LaDa :: load_n_save;
          typedef AtomFreezeMixin  t1;
          lns::xpr::Section main = 
                 (
                   lns::section("Atom")
                     << lns::option( "pos", lns::action=pos,
                                     lns::help="Cartesian position in Anstrom." )
                     << lns::option( "type", lns::action=type,
                                     lns::help="Atomic specie, or any string." )
                     << lns::option( "site", lns::action=site,
                                     lns::help="Atomic site w.r.t. to a lattice." )
                 );
           return _ar & lns::merge(main, *static_cast<t1*>(this));
         }
#   endif
    template<class T> std::ostream& operator<<(std::ostream& _stream, AtomData<T>const & _in)
    {
      return _stream << std::fixed << std::setprecision(5)
                     << std::setw(8) << _in.pos[0] << "  "
                     << std::setw(8) << _in.pos[1] << "  " 
                     << std::setw(8) << _in.pos[2] << "  "
                     << details::print_occupation(_in.type);
    }


  } // namespace Crystal
} // namespace LaDa
  
#endif
