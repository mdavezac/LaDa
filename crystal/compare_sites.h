
#ifndef _LADA_COMPARE_SITES_H_
#define _LADA_COMPARE_SITES_H_

#include "LaDaConfig.h"

#include <set>

#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/mpl/or.hpp>
#include <boost/foreach.hpp>

#include <set>
#include "atom.h"
#include "is_container.h"
#include "utilities.h"


namespace LaDa
{
  namespace crystal 
  {

    namespace details
    {
      //! \brief A functor to compare atomic sites according to type and position.
      //! \details Types and positions can be compared jointly or separately.
      //!          Positions are compared as is without translational
      //!          invariance.
      template<class T_TYPE, class ENABLE=void> struct CompareSites;

      //! Functor to compare sites with container occupations.
      template<class T_TYPE> 
        struct CompareSites<T_TYPE, typename boost::disable_if< is_container<T_TYPE> >::type>
        {
          template <int> struct dummy { dummy(int) {} };
          //! Constructor.
          CompareSites(Atom<T_TYPE> const &_site, types::t_real _tolerance = -1e0 )
          {
            tolerance = _tolerance < 0e0 ? types::tolerance: _tolerance;
            pos = _site.pos;
            type = _site.type;
          }
          //! Copy constructor.
          CompareSites( CompareSites const &_c ): type(_c.type), pos(_c.pos), tolerance(_c.tolerance) {}
          //! Compares string occupation.
          bool operator()( T_TYPE const& _type ) const { return cmp_<T_TYPE>(_type); }
          //! Compares positions.
          bool operator()( math::rVector3d const& _pos ) const { return math::eq(_pos, pos, tolerance); }
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
              cmp_(T const &_type) const { return eq(_type, type, T(tolerance)); }
          math::rVector3d pos;
          types::t_real tolerance;
          T_TYPE type;
        };

      //! Functor to compare sites with scalar occupations.
      template<class T_TYPE> 
        struct CompareSites<T_TYPE, typename boost::enable_if< is_container<T_TYPE> >::type>
        {
          //! Constructor.
          CompareSites(Atom<T_TYPE> const &_site, types::t_real _tolerance = -1e0 )
          {
            tolerance = _tolerance < 0e0 ? types::tolerance: _tolerance;
            pos = _site.pos;
            std::copy(_site.type.begin(), _site.type.end(), std::inserter(type, type.begin()));
          }
          //! Copy constructor.
          CompareSites( CompareSites const &_c ): type(_c.type), pos(_c.pos), tolerance(_c.tolerance) {}
          //! Compares occupations.
          bool operator()( T_TYPE const& _type ) const { return cmp_<T_TYPE>(_type); }
          //! Compares positions.
          bool operator()( math::rVector3d const& _pos ) const { return math::eq(_pos, pos, tolerance); }
          //! Compares positions and site occupation.
          bool operator()( Atom<T_TYPE> const& _site ) const 
            { return operator()(_site.type) and operator()(_site.pos); }

          private:
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
                  if( type.find(_t) != type.end() ) return false;
                return true;
              }
          //! Compares floating point occupations.
          template<class T>
            typename boost::enable_if< boost::is_floating_point<T>, bool>::type
              cmp_( T const _type ) const
              {
                typename std::set<typename T::value_type> :: const_iterator const i_end = type.end();
                foreach(typename T::value_type const &_t, _type)
                {
                  typename std::set<typename T::value_type> :: const_iterator i_first = type.begin();
                  for(; i_first != i_end; ++i_first)
                    if( neq(*i_first, _t, tolerance) ) return false;
                }
                return true;
              }
          math::rVector3d pos;
          types::t_real tolerance;
          std::set<typename T_TYPE::value_type> type;
        };
    }

    //! Returns a functor to compare sites according to positions and/or occupations.
    template<class T_TYPE>
      details::CompareSites<T_TYPE> compare_sites(Atom<T_TYPE> const &_origin, types::t_real _tolerance = -1e0)
        { return details::CompareSites<T_TYPE>(_origin, _tolerance); }


  } // namespace Crystal
} // namespace LaDa
#endif
