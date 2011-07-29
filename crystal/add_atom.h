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

#define BOOST_ASSIGN_MAX_PARAMS 8
#include <boost/assign/list_inserter.hpp>

#include "atom.h"
#include "is_container.h"

namespace LaDa
{
  namespace crystal
  {
    namespace details
    {
      //! Allows easy structure initialization.
      template<class T_TYPE, class ENABLE=void>  struct AddAtom;
      //! Allows easy structure initialization for non-container types.
      template<class T_TYPE>  
        struct AddAtom<T_TYPE,
           typename boost::disable_if< is_container<typename T_TYPE::t_Type> >::type >
        {
          //! type of the species.
          typedef typename T_TYPE::t_Type specie_type;
          //! Atom container.
          std::vector<T_TYPE> &container_;
          //! Constructor.
          AddAtom(std::vector<T_TYPE> &_in) : container_(_in) {}
          //! adder itself.
          void operator()(types::t_real _x, types::t_real _y, types::t_real _z, specie_type const &_type)
          { 
            T_TYPE atom;
            atom.pos[0] = _x; atom.pos[1] = _y; atom.pos[2] = _z; 
            atom.type = _type;
            container_.push_back(_type);
          }
          void operator()(T_TYPE const &_in) { container_.push_back(_in); }
        };
      //! Allows easy structure initialization.
      template<class T_TYPE>  
        struct AddAtom<T_TYPE, 
           typename boost::enable_if< is_container<typename T_TYPE::t_Type> >::type >
        {
          //! type of the species.
          typedef typename T_TYPE::t_Type specie_type;
          //! Atom container.
          std::vector<T_TYPE> &container_;
          //! Constructor.
          AddAtom(std::vector<T_TYPE> &_in) : container_(_in) {}
#         ifdef LADA_TEXT
#           error LADA_TEXT already defined
#         endif
#         ifdef LADA_OPERATOR
#           error LADA_OPERATOR already defined
#         endif
#         define LADA_TEXT(z, n, text) atom.type.push_back(_t ## n);
#         define LADA_OPERATOR(z, n, _)\
            void operator()( types::t_real _x, types::t_real _y, types::t_real _z,      \
                             BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), specie_type const &_t) )  \
            {                                                                           \
              T_TYPE atom;                            \
              atom.pos[0] = _x; atom.pos[1] = _y; atom.pos[2] = _z;                     \
              BOOST_PP_REPEAT_ ## z(BOOST_PP_INC(n), LADA_TEXT, nil)                    \
              container_.push_back(atom);                                               \
            }
          BOOST_PP_REPEAT(5, LADA_OPERATOR, nil)
#         undef LADA_TEXT
#         undef LADA_OPERATOR
          void operator()(T_TYPE const &_in) { container_.push_back(_in); }
        };

      template<class T_DERIVED>
        struct call_add_atom_base
        {
          typedef AddAtom<typename T_DERIVED::value_type> t_Add; 
          boost::assign::list_inserter<t_Add>
            add_atom(typename T_DERIVED::value_type const &_in)
            {
              namespace ba = boost::assign;
              typename T_DERIVED::t_Atoms const &atoms =
                static_cast<T_DERIVED*>(this)->atoms;
              return ba::list_inserter<t_Add>(t_Add(atoms))(_in);
            }
        };

      template<class T_DERIVED, class ENABLE=void> struct call_add_atom;
      template<class T_DERIVED>
        struct call_add_atom<  
                               T_DERIVED, 
                               typename boost::disable_if
                                 <is_container<typename T_DERIVED::t_Type> >::type
                            > :  public call_add_atom_base<T_DERIVED>
        {
          typedef typename T_DERIVED::t_Type specie;
//         typedef AddAtom<typename T_DERIVED::value_type> t_Add; 
          // boost::assign::list_inserter<t_Add>
          void add_atom( types::t_real _x, types::t_real _y, types::t_real _z, specie const &_t)
            {}
          //   namespace ba = boost::assign;
          //   typename T_DERIVED::t_Atoms const &atoms =
          //     static_cast<T_DERIVED*>(this)->atoms;
          //   return ba::list_inserter<t_Add>(t_Add(atoms))(_x, _y, _z,_t);
          // }
        };
      template<class T_DERIVED>
        struct call_add_atom<
                               T_DERIVED, 
                               typename boost::enable_if
                                   <is_container<typename T_DERIVED::value_type::t_Type> >::type
                            > :  public call_add_atom_base<T_DERIVED>
        {
          typedef typename T_DERIVED::value_type::t_Type specie;
          typedef AddAtom<typename T_DERIVED::value_type> t_Add; 
#         define LADA_TEXT(z, n, text) atom.type.push_back(_t ## n);
#         define LADA_OPERATOR(z, n, _)                                               \
            boost::assign::list_inserter<t_Add> \
              add_atom( types::t_real _x, types::t_real _y, types::t_real _z,         \
                        BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), specie const &_t) )     \
              {                                                                       \
                namespace ba = boost::assign;                                         \
                typename T_DERIVED::t_Atoms const &atoms =                    \
                  static_cast<T_DERIVED*>(this)->atoms;                       \
                return ba::list_inserter<t_Add>(t_Add(atoms))   \
                         ( _x, _y, _z,                                                \
                           BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), _t) );               \
              }
          BOOST_PP_REPEAT(5, LADA_OPERATOR, nil)
#         undef LADA_TEXT
#         undef LADA_OPERATOR
        };

      template<class T_DERIVED> struct call_add_atom2_base
      {
        typedef AddAtom<typename T_DERIVED::value_type> t_Add; 
        boost::assign::list_inserter<t_Add>
          add_atom(typename T_DERIVED::value_type const &_in)
            { return static_cast<T_DERIVED*>(this)->impl_->add_atom(_in); }
      };
      template<class T_DERIVED, class ENABLE=void> struct call_add_atom2;
      template<class T_DERIVED>
        struct call_add_atom2< 
                                T_DERIVED,
                                typename boost::disable_if
                                   < is_container<typename T_DERIVED::value_type::t_Type> >::type
                             > : public call_add_atom2_base<T_DERIVED>
        {
          //! type of the species.
          typedef typename T_DERIVED::value_type::t_Type specie;
          typedef AddAtom<typename T_DERIVED::value_type> t_Add; 
          boost::assign::list_inserter<t_Add>
            add_atom( types::t_real _x, types::t_real _y, types::t_real _z, specie const &_t)
              { return static_cast<T_DERIVED*>(this)->impl_->add_atom(_x, _y, _z, _t); }
        };
      template<class T_DERIVED>
        struct call_add_atom2<
                                T_DERIVED, 
                                typename boost::enable_if
                                   < is_container<typename T_DERIVED::value_type::t_Type> >::type
                             > : public call_add_atom2_base<T_DERIVED>
        {
          //! type of the species.
          typedef typename T_DERIVED::value_type::t_Type specie;
          typedef AddAtom<typename T_DERIVED::value_type> t_Add; 
#         define LADA_TEXT(z, n, text) atom.type.push_back(_t ## n);
#         define LADA_OPERATOR(z, n, _)                                                 \
            boost::assign::list_inserter<t_Add> \
              add_atom( types::t_real _x, types::t_real _y, types::t_real _z,           \
                        BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), specie const &_t) )      \
              {                                                                         \
                return static_cast<T_DERIVED*>(this)->impl_                     \
                   ->add_atom(_x, _y, _z, BOOST_PP_ENUM_PARAMS(BOOST_PP_INC(n), _t));   \
              }
          BOOST_PP_REPEAT(5, LADA_OPERATOR, nil)
#         undef LADA_TEXT
#         undef LADA_OPERATOR
        };
       
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

