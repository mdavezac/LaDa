#ifndef LADA_CRYSTAL_ADDATOM_H
#define LADA_CRYSTAL_ADDATOM_H

#include "LaDaConfig.h"

#include <iterator>

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
#     ifdef LADA_BASIC
#       error LADA_BASIC already defined
#     endif
#     ifdef LADA_TEXT
#       error LADA_TEXT already defined
#     endif
#     ifdef LADA_TEXT1
#       error LADA_TEXT1 already defined
#     endif
#     ifdef LADA_OPERATOR
#       error LADA_OPERATOR already defined
#     endif
      template<class T_TYPE, class Enable=void> struct AddAtom;

#     define LADA_BASIC                                                                             \
        typedef T_TYPE t_Type;                                                                      \
        typedef typename traits::StructureData<t_Type>::t_Atom t_Atom;                              \
        typedef typename traits::StructureData<t_Type>::t_Atoms t_Atoms;                            \
        typedef typename std::back_insert_iterator<t_Atoms> t_Inserter;                             \
        AddAtom(t_Inserter const &_inserter) : inserter_(_inserter) {}                              \
        AddAtom(AddAtom const &_c) : inserter_(_c.inserter_) {}                                     \
        AddAtom operator()(types::t_real _x, types::t_real _y, types::t_real _z, T_TYPE const &_t)  \
          { return operator()(math::rVector3d(_x, _y, _z), _t); }                                   \
        template<class TD>                                                                          \
          AddAtom operator()(Eigen::DenseBase<TD> const &_pos, T_TYPE const &_t)                    \
            { return operator()(t_Atom(_pos, _t)); }                                                \
        AddAtom operator()(t_Atom const& _in) { *inserter_ = _in; ++inserter_; return *this; }      \
        private:                                                                                    \
         t_Inserter inserter_;
       

      //! Actually performs atom addition.
      template<class T_TYPE> struct AddAtom
        <
          T_TYPE,
          typename boost::enable_if< details::is_scalar<T_TYPE> >::type
        > 
        {
          LADA_BASIC;
        };

      //! Actually performs atom addition, specialize for iterables.
      template<class T_TYPE> struct AddAtom
        <
          T_TYPE,
          typename boost::enable_if< details::is_iterable<T_TYPE> >::type
        >
      {
        LADA_BASIC;
        public:
        typedef typename is_container<t_Type>::type cont_or_set;
        typedef typename t_Type::const_reference specie_cref;
#       define LADA_OPERATOR(z, n, _)                                                        \
          AddAtom operator()( types::t_real _x, types::t_real _y, types::t_real _z,          \
                              BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), specie_cref _t) )        \
            {                                                                                \
              t_Type type;                                                                   \
              fill(type, cont_or_set(), BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), _t));          \
              return operator()(math::rVector3d(_x, _y, _z), type);                          \
            }                                                                                \
          template<class TD>                                                                 \
            AddAtom operator()( Eigen::DenseBase<TD> const &_pos,                            \
                                BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), specie_cref _t) )      \
              {                                                                              \
                t_Type type;                                                                 \
                fill(type, cont_or_set(), BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), _t));        \
                return operator()(_pos, type);                                               \
              }
        BOOST_PP_REPEAT(5, LADA_OPERATOR, nil)
#       undef LADA_OPERATOR
        private:
#       define LADA_TEXT(z, n, text) _in.push_back(_t ## n);
#       define LADA_TEXT1(z, n, text) _in.insert(_t ## n);
#       define LADA_OPERATOR(z, n, _)                                                      \
          void fill( t_Type& _in, boost::mpl::bool_<true>,                                 \
                     BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), specie_cref _t) )               \
          {                                                                                \
            BOOST_PP_REPEAT_ ## z(BOOST_PP_INC(n), LADA_TEXT, nil)                         \
          }                                                                                \
          void fill( t_Type& _in, boost::mpl::bool_<false>,                                \
                     BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), specie_cref _t) )               \
          {                                                                                \
            BOOST_PP_REPEAT_ ## z(BOOST_PP_INC(n), LADA_TEXT1, nil)                        \
          }
        BOOST_PP_REPEAT(5, LADA_OPERATOR, nil)
#       undef LADA_OPERATOR
#       undef LADA_TEXT
#       undef LADA_TEXT1
      };

      //! Base of the mixin for adding atoms to a structure.
      template<template <class> class T_DERIVED, class T_TYPE> 
        struct AddAtomMixinBase
        {
          typedef typename traits::StructureData<T_TYPE>::t_Atom t_Atom;
          typedef typename traits::StructureData<T_TYPE>::t_Atoms t_Atoms;
          AddAtom<T_TYPE> add_atom( types::t_real _x, types::t_real _y,
                                    types::t_real _z, T_TYPE const &_t) 
          { return add_atom(math::rVector3d(_x, _y, _z), _t); }
          template<class TD>                                                                    
            AddAtom<T_TYPE> add_atom(Eigen::DenseBase<TD> const &_pos, T_TYPE const &_t)      
              { return add_atom(t_Atom(_pos, _t)); }                                          
          AddAtom<T_TYPE> add_atom(t_Atom const& _in)
            { return AddAtom<T_TYPE>(this->back_inserter())(_in); }
          protected:
          std::back_insert_iterator<t_Atoms> back_inserter() 
            { return std::back_inserter(static_cast<T_DERIVED<T_TYPE>*>(this)->atoms()); }
        };
      
      //! \brief Class to easily add an atom to a StructureData instance.
      //! \details We want three different points.
      //!            - We should be able to distinguish between structures
      //!              where the type is a scalar and those where it is a container (or a set).
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
      template<template <class> class T_DERIVED, class T_TYPE, class ENABLE=void> struct AddAtomMixin;
      //! Specializatio for scalar types.
      template<template <class> class T_DERIVED, class T_TYPE> 
        struct AddAtomMixin<  
                             T_DERIVED, T_TYPE,
                             typename boost::enable_if< is_scalar<T_TYPE> >::type
                           > : public AddAtomMixinBase<T_DERIVED, T_TYPE> {};
      //! Specializatio for container and set types.
      template<template <class> class T_DERIVED, class T_TYPE> 
        struct AddAtomMixin<  
                             T_DERIVED, T_TYPE,
                             typename boost::enable_if< is_iterable<T_TYPE> >::type
                           > : public AddAtomMixinBase<T_DERIVED, T_TYPE> 
        {
          typedef typename traits::StructureData<T_TYPE>::t_Atom t_Atom;
          typedef typename T_TYPE::value_type specie;
#         define LADA_OPERATOR(z, n, _)                                                    \
            AddAtom<T_TYPE>                                                                \
              add_atom( types::t_real _x, types::t_real _y, types::t_real _z,              \
                        BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), specie const &_t) )          \
              { return add_atom( math::rVector3d(_x, _y, _z),                              \
                                 BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), _t) ); }            \
            template<class TD>                                                             \
              AddAtom<T_TYPE>                                                              \
                add_atom( Eigen::DenseBase<TD> const &_pos,                                \
                          BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), specie const &_t) )        \
                { return AddAtom<T_TYPE>(this->back_inserter())                            \
                     (_pos, BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), _t) ); }                  
            BOOST_PP_REPEAT(5, LADA_OPERATOR, nil)
#         undef LADA_OPERATOR
        };
#     undef LADA_BASIC

    } //namespace details
  } // namespace Crystal
} // namespace LaDa
  
#endif

