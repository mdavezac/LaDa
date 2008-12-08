//
//  Version: $Id: taboos.h 856 2008-11-14 17:00:23Z davezac $
//
#ifndef _DARWIN_TABOO_H_
#define _DARWIN_TABOO_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <print/xmg.h>

namespace LaDa
{
  namespace GA
  {
    namespace Taboo
    {
      template< class T_INDIVIDUAL, class T_CONTAINER = std::vector< T_INDIVIDUAL > >
        class History 
        {
          public:
            //! Type of the population container.
            typedef T_CONTAINER t_Container;
            //! Type of the individual
            typedef T_INDIVIDUAL t_Individual;
            //! Mimics container stuff.
            typedef t_Individual value_type;
            //! Type of the iterators.
            typedef typename t_Container :: const_iterator const_iterator;


            //! Constructor.
            History() {} 
            //! Copy Constructor.
            History( const History& _c )
            {
              if( _c.is_on() ) container_ = _c.container_; 
              else off();
            }
            //! Destructor.
            virtual ~History() {}

            //! The functor. 
            bool operator()( const t_Individual& _indiv ) const
            {
              if( not is_on() ) return false;
              return container_->end() == std::find( container_->begin(), container_->end(), _indiv ); 
            }
            
            //! If \a _indiv exists in History::container_, copies it.
            bool clone(t_Individual &_indiv) const;
            //! Adds an individual to the taboo list, if it does not exist.
            bool push_back( t_Individual &_indiv )
            {
              if( not is_on() ) return false;
              if( not (*this)( _indiv ) ) container_->push_back( _indiv ); 
            }

            //! Turn history on.
            void on() { if( not is_on() ) container_.reset( new t_Container ); }
            //! Turn history off.
            void off() { if( is_on() ) container_.reset(); }
            //! returns true if history is on.
            bool is_on() const { return bool(container_); }
            //! returns true if history is off.
            bool is_off() const { return not bool(container_); }
            //! Number of individuals in history.
            size_t size() const
              { return is_off() ? 0: container_->size(); }
            //! Iterator to container. Undefined if history is off.
            typename t_Container :: const_iterator begin() const 
              { return  container_->begin(); }
            //! Iterator to container. Undefined if history is off.
            typename t_Container :: const_iterator end() const 
              { return  container_->end(); }


          protected:
            //! A reference to the islands.
            boost::shared_ptr< t_Container > container_;
        };

      template< class T_INDIVIDUAL, class T_CONTAINER >
        bool History<T_INDIVIDUAL, T_CONTAINER> :: clone( t_Individual& _indiv ) const
        {
          if( not is_on() ) return false;
          typedef typename t_Container :: const_iterator t_cit;
          t_cit i_found = std::find( container_->begin(), container_->end(), _indiv );
          if( i_found == container_->end() ) return false;
          _indiv = *i_found;
          return true;
        }

      namespace Factory
      {
        //! Factory for creating a container of taboos.
        template< class T_FACTORY, class T_INDIVIDUAL >
          void container( T_FACTORY &_factory,
                          boost::function<bool( const T_INDIVIDUAL& )>& _function,
                          const TiXmlElement &_node,
                          Taboo::History<T_INDIVIDUAL>& _history )
          {
            Print::xmg << Print::Xmg::comment << "History." << Print::endl;
            _history.on();
            _function = _history; 
          }
      }

    } // namespace Taboo
  } // namespace GA
} // namespace LaDa
#endif
