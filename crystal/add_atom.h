#ifndef LADA_CRYSTAL_ADDATOM_H
#define LADA_CRYSTAL_ADDATOM_H

#include "LaDaConfig.h"

#include <boost/serialization/base_object.hpp>
#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/enum_params.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/next.hpp>
#include <boost/mpl/int.hpp>

#include "is_container.h"
#include "traits.h"

namespace LaDa
{
  namespace crystal
  {
    namespace details
    {
#     ifdef LADA_TEXT
#       error LADA_TEXT already defined
#     endif
#     ifdef LADA_OPERATOR
#       error LADA_OPERATOR already defined
#     endif
      //! \brief Class to easily add an atom to a StructureData instance.
      //! \details We want three different points.
      //!            - We should be able to distinguish between structures
      //!              where the type is a scalar and those where it is a container.
      //!            - When considering scalar types, the structure should only
      //!              contain add_atom functions with 4 arguments (3
      //!              positions, 1 type). When considering a container type,
      //!              the structure should contain multiple overloaded add_atoms functions.
      //!            - add_atom should be member functions.
      //!            .
      //!          To deal with the first two points, eg scalar vs container
      //!          issues, we use boost::enable_if to instantiate only the
      //!          correct version of add_atom. To deal the third point, we
      //!          incorporate the functions using the curriously recurring
      //!          template idiom.
      template<template <class> class T_DERIVED, class T_TYPE, class ENABLE=void> struct call_add_atom;
      //! Specializatio for scalar types and StructureData.
      template<template <class> class T_DERIVED, class T_TYPE> 
        struct call_add_atom<  
                               T_DERIVED, T_TYPE,
                               typename boost::disable_if< is_container<T_TYPE> >::type
                            > 
        {
          typedef typename traits::StructureData<T_TYPE>::t_Atom t_Atom;
          call_add_atom& add_atom(types::t_real _x, types::t_real _y, types::t_real _z, T_TYPE const &_t)
          {
            t_Atom atom;
            atom.pos[0] = _x; atom.pos[1] = _y; atom.pos[2] = _z; 
            atom.type = _t;
            static_cast<T_DERIVED<T_TYPE>*>(this)->atoms.push_back(atom);
            return *this;
          }
          call_add_atom& add_atom(t_Atom const &_in)
          {
            static_cast<T_DERIVED<T_TYPE>*>(this)->atoms.push_back(_in);
            return *this;
          }
        };
      //! Specializatio for container types and StructureData.
      template<template<class> class T_DERIVED, class T_TYPE>
        struct call_add_atom<
                               T_DERIVED, T_TYPE,
                               typename boost::enable_if< is_container<T_TYPE> >::type
                            >
        {
          typedef typename traits::StructureData<T_TYPE>::t_Atom t_Atom;
          typedef typename T_TYPE::value_type specie;
#         define LADA_TEXT(z, n, text) atom.type.push_back(_t ## n);
#         define LADA_OPERATOR(z, n, _)                                               \
            call_add_atom&                                                            \
              add_atom( types::t_real _x, types::t_real _y, types::t_real _z,         \
                        BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), specie const &_t) )     \
              {                                                                       \
                t_Atom atom;                                                          \
                atom.pos[0] = _x; atom.pos[1] = _y; atom.pos[2] = _z;                 \
                BOOST_PP_REPEAT_ ## z(BOOST_PP_INC(n), LADA_TEXT, nil)                \
                static_cast<T_DERIVED<T_TYPE>*>(this)->atoms.push_back(atom);         \
                return *this;                                                         \
              }
          BOOST_PP_REPEAT(5, LADA_OPERATOR, nil)
#         undef LADA_TEXT
#         undef LADA_OPERATOR
          call_add_atom& add_atom(t_Atom const &_in)
          {
            static_cast<T_DERIVED<T_TYPE>*>(this)->atoms.push_back(_in);
            return *this;
          }
        };

      //! Replicates class call_add_atom for TemplateStructure.
      template<template <class> class T_DERIVED, class T_TYPE, class ENABLE=void> struct call_add_atom2;
      //! Specializatio for scalar types and TemplateStructure.
      template<template <class> class T_DERIVED, class T_TYPE> 
        struct call_add_atom2<  
                                T_DERIVED, T_TYPE,
                                typename boost::disable_if< is_container<T_TYPE> >::type
                             > 
        {
          typedef typename traits::StructureData<T_TYPE>::t_Atom t_Atom;
          call_add_atom2& add_atom( types::t_real _x, types::t_real _y, types::t_real _z,
                                    T_TYPE const &_t )
          {
            static_cast<T_DERIVED<T_TYPE>*>(this)->impl_->add_atom(_x, _y, _z, _t);
            return *this;
          }
          call_add_atom2& add_atom(t_Atom const &_in)
          {
            static_cast<T_DERIVED<T_TYPE>*>(this)->impl_->add_atom(_in);
            return *this;
          }
        };
      //! Specializatio for container types and TemplateStructure.
      template<template<class> class T_DERIVED, class T_TYPE>
        struct call_add_atom2<
                                T_DERIVED, T_TYPE,
                                typename boost::enable_if< is_container<T_TYPE> >::type
                             >
        {
          typedef typename traits::StructureData<T_TYPE>::t_Atom t_Atom;
          typedef call_add_atom2<T_DERIVED, T_TYPE> &t_Return;
          typedef typename T_TYPE::value_type specie;
#         define LADA_TEXT(z, n, text) atom.type.push_back(_t ## n);
#         define LADA_OPERATOR(z, n, _)                                                  \
            t_Return add_atom( types::t_real _x, types::t_real _y, types::t_real _z,     \
                               BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), specie const &_t) ) \
              {                                                                          \
                static_cast<T_DERIVED<T_TYPE>*>(this)->                                  \
                  impl_->add_atom( _x, _y, _z,                                           \
                                   BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), _t) );          \
                return *this;                                                            \
              }
          BOOST_PP_REPEAT(5, LADA_OPERATOR, nil)
#         undef LADA_TEXT
#         undef LADA_OPERATOR
          t_Return add_atom(t_Atom const &_in)
          {
            static_cast<T_DERIVED<T_TYPE>*>(this)->impl_->add_atom(_in);
            return *this;
          }
        };

      //! Easily sets the cell parameters of a StructureData and TemplateData.
      template<class D = boost::mpl::int_<0> >
        class SetCell
        {
          public:
            SetCell(math::rMatrix3d &_cell) : cell_(_cell) {}
            SetCell<typename boost::mpl::next<D>::type>
              operator()(types::t_real _x, types::t_real _y, types::t_real _z)
              {
                cell_(D::value, 0) = _x;
                cell_(D::value, 1) = _y;
                cell_(D::value, 2) = _z;
                return SetCell<typename boost::mpl::next<D>::type>(cell_);
              }
            SetCell<typename boost::mpl::next<D>::type>
              operator()(math::rVector3d const &_pos)
              {
                cell_(D::value, 0) = _pos[0];
                cell_(D::value, 1) = _pos[1];
                cell_(D::value, 2) = _pos[2];
                return SetCell<typename boost::mpl::next<D>::type>(cell_);
              }
          private:
            math::rMatrix3d &cell_;
        };
      //! \brief Easily sets the cell parameters of a StructureData and TemplateData.
      //! \details This one returns void, so that only three instances of
      //!          SetCell can be called.
      template<> class SetCell< boost::mpl::int_<2> >
      {
        public:
          SetCell(math::rMatrix3d &_cell) : cell_(_cell) {}
          void operator()(types::t_real _x, types::t_real _y, types::t_real _z)
          {
            cell_(boost::mpl::int_<2>::value, 0) = _x;
            cell_(boost::mpl::int_<2>::value, 1) = _y;
            cell_(boost::mpl::int_<2>::value, 2) = _z;
          }
          void operator()(math::rVector3d const &_pos)
          {
            cell_(boost::mpl::int_<2>::value, 0) = _pos[0];
            cell_(boost::mpl::int_<2>::value, 1) = _pos[1];
            cell_(boost::mpl::int_<2>::value, 2) = _pos[2];
          }
        private:
          math::rMatrix3d &cell_;
      };
    }

  } // namespace Crystal
} // namespace LaDa
  
#endif

