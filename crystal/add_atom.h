#ifndef LADA_CRYSTAL_ATOM_H
#define LADA_CRYSTAL_ATOM_H

#include "LaDaConfig.h"

#include <boost/serialization/base_object.hpp>
#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/enum_params.hpp>
#include <boost/utility/enable_if.hpp>

#define BOOST_ASSIGN_MAX_PARAMS 8
#include <boost/assign/list_inserter.hpp>

#include "atom.h"

namespace LaDa
{
  namespace crystal
  {
    namespace details
    {
      template<class T_CONTAINER>
        struct is_container : public boost::mpl::false_ {};
      template< template<class A, class B> class T_CONTAINER, class TYPE, class ALLOC>
        struct is_container<T_CONTAINER<TYPE, ALLOC> >: public boost::mpl::true_ {};
      //! Allows easy structure initialization.
      template<class T_TYPE, class ENABLE>  struct add_atom;
      //! Allows easy structure initialization for non-container types.
      template<class T_TYPE>  
        struct add_atom<T_TYPE, boost::disable_if< typename is_container<T_TYPE>::type >
        {
          //! Atom container.
          std::vector<T_TYPE> &container_;
          //! Constructor.
          add_atom(std::vector<T_TYPE> &_in) : container_(_in) {}
          //! adder itself.
          void operator()( types::t_real _x, types::t_real _y, types::t_real _z, T_TYPE const &_type)
          { 
            Atom<T_TYPE> atom;
            atom.pos[0] = _x; atom.pos[1] = _y; atom.pos[2] = _z; 
            atom.type = _type;
            container_.push_back(_type);
          }
          void operator()(Atom<TYPE> const &_in) { atoms.push_back(_in); }
        };
      //! Allows easy structure initialization.
      template<class T_TYPE>  
        struct add_atom<T_TYPE, boost::enable_if< typename is_container<T_TYPE>::type >
        {
          //! Atom container.
          std::vector<T_TYPE> &container_;
          //! Constructor.
          add_atom(std::vector<T_TYPE> &_in) : container_(_in) {}
#         ifdef LADA_TEXT
#           error LADA_TEXT already defined
#         endif
#         ifdef LADA_OPERATOR
#           error LADA_OPERATOR already defined
#         endif
#         define LADA_TEXT(z, n, text) atom.type.push_back(_t ## n);
#         define LADA_OPERATOR(z, n, _)\
            void operator()( types::t_real _x, types::t_real _y, types::t_real _z,      \
                             BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), t_Type const &_t) )  \
            {                                                                           \
              Atom<T_TYPE> atom;                                                        \
              atom.pos[0] = _x; atom.pos[1] = _y; atom.pos[2] = _z;                     \
              BOOST_PP_REPEAT_ ## z(BOOST_PP_INC(n), LADA_TEXT, nil)                    \
              container_.push_back(atom);                                               \
            }
          BOOST_PP_REPEAT(5, LADA_OPERATOR, nil)
#         undef LADA_TEXT
#         undef LADA_OPERATOR
          void operator()(Atom<TYPE> const &_in) { atoms.push_back(_in); }
        };

      template<template<class> class T_DERIVED, class T_TYPE, class ENABLE> struct call_add_atom
      template<template<class> class T_DERIVED, class T_TYPE> 
        struct call_add_atom<T_DERIVED, T_TYPE, typename boost::disable_if< is_container<T_TYPE>::type >
        {
            boost::assign::list_inserter<details::add_atom>
              add_atom( types::t_real _x, types::t_real _y, types::t_real _z, T_TYPE const &_t)
              {
                namespace ba = boost::assign;
                typename T_DERIVED<T_TYPE>::t_Atoms const &atoms =
                  static_cast<T_DERIVED<T_TYPE>*>(this)->atoms;
                return ba::list_inserter(details::add_atom(atoms))(_x, _y, _z,_t);
              }
            boost::assign::list_inserter<details::add_atom>
              add_atom(Atom<T_TYPE> const &_in)
              {
                namespace ba = boost::assign;
                typename T_DERIVED<T_TYPE>::t_Atoms const &atoms =
                  static_cast<T_DERIVED<T_TYPE>*>(this)->atoms;
                return ba::list_inserter(details::add_atom(atoms))(_in);
              }
        };
      template<template<class> class T_DERIVED, class T_TYPE> 
        struct call_add_atom<T_DERIVED, T_TYPE, typename boost::enable_if< is_container<T_TYPE>::type >
        {
#         define LADA_TEXT(z, n, text) atom.type.push_back(_t ## n);
#         define LADA_OPERATOR(z, n, _)\
            boost::assign::list_inserter<details::add_atom>                           \
              add_atom( types::t_real _x, types::t_real _y, types::t_real _z,         \
                        BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), t_Type const &_t) )     \
              {                                                                       \
                namespace ba = boost::assign;                                         \
                typename T_DERIVED<T_TYPE>::t_Atoms const &atoms =                    \
                  static_cast<T_DERIVED<T_TYPE>*>(this)->atoms;                       \
                return ba::list_inserter(details::add_atom(atoms))                    \
                         ( _x, _y, _z,                                                \
                           BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), _t) );               \
              }
          BOOST_PP_REPEAT(5, LADA_OPERATOR, nil)
#         undef LADA_TEXT
#         undef LADA_OPERATOR
            boost::assign::list_inserter<details::add_atom>
              add_atom(Atom<T_TYPE> const &_in)
              {
                namespace ba = boost::assign;
                typename T_DERIVED<T_TYPE>::t_Atoms const &atoms =
                  static_cast<T_DERIVED<T_TYPE>*>(this)->atoms;
                return ba::list_inserter(details::add_atom(atoms))(_in);
              }
        };

      template<template<class> class T_DERIVED, class T_TYPE, class ENABLE> struct call_add_atom2
      template<template<class> class T_DERIVED, class T_TYPE> 
        struct call_add_atom2<T_DERIVED, T_TYPE, typename boost::disable_if< is_container<T_TYPE>::type >
        {
          boost::assign::list_inserter<details::add_atom>
            add_atom( types::t_real _x, types::t_real _y, types::t_real _z, T_TYPE const &_t)
              { return static_cast< T_DERIVED<T_TYPE> >(*this)->impl_->add_atom(_x, _y, _z, _t); }
          boost::assign::list_inserter<details::add_atom>
            add_atom(Atom<T_TYPE> const &_in)
              { return static_cast< T_DERIVED<T_TYPE> >(*this)->impl_->add_atom(_in); }
        };
      template<template<class> class T_DERIVED, class T_TYPE> 
        struct call_add_atom2<T_DERIVED, T_TYPE, typename boost::enable_if< is_container<T_TYPE>::type >
        {
#         define LADA_TEXT(z, n, text) atom.type.push_back(_t ## n);
#         define LADA_OPERATOR(z, n, _)\
            boost::assign::list_inserter<details::add_atom>                             \
              add_atom( types::t_real _x, types::t_real _y, types::t_real _z,           \
                        BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), t_Type const &_t) )       \
              {                                                                         \
                return static_cast< T_DERIVED<T_TYPE> >(*this)->impl_                   \
                   ->add_atom(_x, _y, _z, BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), _t));   \
              }
          BOOST_PP_REPEAT(5, LADA_OPERATOR, nil)
#         undef LADA_TEXT
#         undef LADA_OPERATOR
          boost::assign::list_inserter<details::add_atom>
            add_atom(Atom<T_TYPE> const &_in)
              { return static_cast< T_DERIVED<T_TYPE> >(*this)->impl_->add_atom(_in); }
        };
       
    }

  } // namespace Crystal
} // namespace LaDa
  
#endif

