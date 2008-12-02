//
//  Version: $Id: taboos.h 856 2008-11-14 17:00:23Z davezac $
//
#ifndef _LADA_GA_TABOO_POPULATIONS_H_
#define _LADA_GA_TABOO_POPULATIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <print/xmg.h>

namespace LaDa
{
  namespace GA
  {
    namespace Taboo
    {
      //! \brief An offspring taboo.
      //! \details Checks in range [0, n-1] where n is the length of \a _population.
      template< class T_POPULATION > 
        bool offspring( const typename T_POPULATION :: value_type& _indiv,
                        const T_POPULATION& _population )
        {
          typename T_POPULATION :: const_iterator i_begin = _population.begin();
          typename T_POPULATION :: const_iterator i_end = _population.end();
          if( i_begin == i_end ) return false;
          --i_end;
          if( i_begin == i_end ) return false;

          return i_end != std::find( i_begin, i_end, _indiv ); 
        }

      //! A population taboo.
      template< class T_POPULATION > 
        bool population( const typename T_POPULATION :: value_type& _indiv,
                         const T_POPULATION& _pop )
          { return _pop.end() != std::find( _pop.begin(), _pop.end(), _indiv ); }

      //! An islands taboo.
      template< class T_ISLANDS > 
        bool islands( const typename T_ISLANDS :: value_type :: value_type& _indiv, 
                      const T_ISLANDS& _islands )
        {
          foreach( const typename T_ISLANDS :: value_type& pop, _islands )
            if( population( _indiv, pop ) ) return true;
          return false;
        }

      //! Factory for taboo operators.
      namespace Factory
      {
        template< class T_CONTAINER >
          void offspring( boost::function<bool( const typename T_CONTAINER::value_type& )>&
                            _function,
                          const TiXmlElement &,
                          const T_CONTAINER& _container )
          {
            _function = boost::bind
                        (
                          &::LaDa::GA::Taboo::offspring<T_CONTAINER>,
                          _1, boost::cref(_container) 
                        );
            Print::xmg << Print::Xmg::comment << "Offspring Taboo" << Print::endl;
          }
        template< class T_CONTAINER >
          void population( boost::function<bool( const typename T_CONTAINER::value_type& )>&
                            _function,
                           const TiXmlElement &,
                           const T_CONTAINER& _container )
          {
            _function = boost::bind
                        (
                          &Taboo::population<T_CONTAINER>,
                          _1, boost::cref(_container) 
                        );
            Print::xmg << Print::Xmg::comment << "Population Taboo" << Print::endl;
          }
        template< class T_CONTAINER >
          void islands( boost::function<bool( const typename T_CONTAINER::value_type::value_type& )>&
                              _function,
                        const TiXmlElement &,
                        const T_CONTAINER& _container )
          {
            _function = boost::bind
                        (
                          &Taboo::islands<T_CONTAINER>,
                          _1, boost::cref(_container) 
                        );
            Print::xmg << Print::Xmg::comment << "Populations Taboo" << Print::endl;
          }
      }
        
    } // namespace Taboo
  } // namespace GA
} // namespace LaDa
#endif
