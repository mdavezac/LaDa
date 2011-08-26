#ifndef LADA_COMPARE_SITES_H
#define LADA_COMPARE_SITES_H

#include "LaDaConfig.h"

#include <algorithm>

#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/mpl/or.hpp>
#include <boost/foreach.hpp>

#include "atom.h"
#include "is_container.h"
#include "utilities.h"


namespace LaDa
{
  namespace crystal 
  {
    //! \brief A functor to compare atomic sites according to type and position.
    //! \details Types and positions can be compared jointly or separately.
    //!          Positions are compared as is without translational
    //!          invariance.
    template<class T_TYPE, class ENABLE=void> struct CompareSites;

    //! Functor to compare sites with container occupations.
    template<class T_TYPE> 
      struct CompareSites<T_TYPE, typename boost::disable_if<details::is_container<T_TYPE> >::type>
      {
        template <int> struct dummy { dummy(int) {} };
        //! Constructor.
        CompareSites(Atom<T_TYPE> const &_site, types::t_real _tolerance = types::tolerance )
        {
          tolerance = _tolerance;
          pos = _site.pos;
          type = _site.type;
        }
        template<class T_DERIVED>
          CompareSites( Eigen::DenseBase<T_DERIVED> const &_pos, 
                        T_TYPE const &_type, types::t_real _tolerance = types::tolerance )
          {
            tolerance = _tolerance;
            pos = _pos;
            type = _type;
          }
        //! Copy constructor.
        CompareSites( CompareSites const &_c ): type(_c.type), pos(_c.pos), tolerance(_c.tolerance) {}
        //! Compares string occupation.
        bool operator()( T_TYPE const& _type ) const { return cmp_<T_TYPE>(_type); }
        //! Compares positions.
        template<class T_DERIVED>
          bool operator()(Eigen::DenseBase<T_DERIVED> const& _pos) const 
            { return math::eq(_pos, pos, tolerance); }
        //! Compares positions and site occupation.
        bool operator()( Atom<T_TYPE> const& _site ) const 
          { return operator()(_site.type) and operator()(_site.pos); }

        private:
        //! Compares string occupation.
        template<class T>
          typename boost::enable_if
            < 
              typename boost::mpl::or_
                < typename boost::is_convertible<T, std::string> :: type,
                  typename boost::is_integral<T> > :: type,
              bool
            >::type cmp_(T const &_type) const { return _type == type; }
        //! Compares floating point occupations.
        template<class T>
          typename boost::enable_if<boost::is_floating_point<T>, bool>::type
            cmp_(T const &_type) const { return math::eq(_type, type, T(tolerance)); }
        math::rVector3d pos;
        types::t_real tolerance;
        T_TYPE type;
      };

    //! Functor to compare sites with scalar occupations.
    template<class T_TYPE> 
      struct CompareSites<T_TYPE, typename boost::enable_if<details::is_container<T_TYPE> >::type>
      {
        //! Constructor.
        CompareSites(Atom<T_TYPE> const &_site, types::t_real _tolerance = types::tolerance )
        {
          tolerance = _tolerance;
          pos = _site.pos;
          std::copy(_site.type.begin(), _site.type.end(), std::back_inserter(type));
        }
        template<class T_DERIVED>
          CompareSites( Eigen::DenseBase<T_DERIVED> const &_pos, 
                        T_TYPE const &_type, types::t_real _tolerance = types::tolerance )
          {
            tolerance = _tolerance;
            pos = _pos;
            std::copy(_type.begin(), _type.end(), std::back_inserter(type));
          }
       
        //! Copy constructor.
        CompareSites( CompareSites const &_c ): type(_c.type), pos(_c.pos), tolerance(_c.tolerance) {}
        //! Compares occupations.
        bool operator()( T_TYPE const& _type ) const { return cmp_<T_TYPE>(_type); }
        //! Compares positions.
        template<class T_DERIVED>
          bool operator()(Eigen::DenseBase<T_DERIVED> const& _pos) const 
            { return math::eq(_pos, pos, tolerance); }
        //! Compares positions and site occupation.
        //! Compares positions and site occupation.
        bool operator()( Atom<T_TYPE> const& _site ) const 
          { return operator()(_site.type) and operator()(_site.pos); }

        private:
        template<class T>
          struct Cmp_
          {
            Cmp_(typename T::value_type const &_a, typename T::value_type _tol) : a(_a), tolerance(_tol) {};
            bool operator()(typename T::value_type const &_b) const
            { return math::eq(a, _b, tolerance); }
            typename T::value_type const & a;
            typename T::value_type tolerance;
          };
        //! Compares string or integer occupations.
        template<class T>
          typename boost::enable_if
            < 
              typename boost::mpl::or_
                < boost::is_convertible<typename T::value_type, std::string>,
                  boost::is_integral<typename T::value_type> > :: type,
              bool
            >::type cmp_( T const &_type ) const
            { 
              foreach(typename T::value_type const &_t, _type)
                if( std::find(type.begin(), type.end(), _t) == type.end() ) return false;
              foreach(typename T_TYPE::value_type const &_t, type)
                if( std::find(_type.begin(), _type.end(), _t) == _type.end() ) return false;
              return true;
            }
        //! Compares floating point occupations.
        template<class T>
          typename boost::enable_if< boost::is_floating_point<typename T::value_type>, bool>::type
            cmp_( T const _type ) const
            {
              foreach(typename T::value_type const &_t, _type)
                if( std::find_if(type.begin(), type.end(), Cmp_<T>(_t, tolerance)) == type.end() ) return false;
              foreach(typename T_TYPE::value_type const &_t, type)
                if( std::find_if(_type.begin(), _type.end(), Cmp_<T>(_t, tolerance)) == _type.end() ) return false;
              return true;
            }
        math::rVector3d pos;
        types::t_real tolerance;
        T_TYPE type;
      };
    

    //! Returns a functor to compare sites according to positions and/or occupations.
    template<class T_TYPE>
      CompareSites<T_TYPE> compare_sites(Atom<T_TYPE> const &_origin, types::t_real _tolerance = types::tolerance)
        { return CompareSites<T_TYPE>(_origin, _tolerance); }


  } // namespace Crystal
} // namespace LaDa
#endif
