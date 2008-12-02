//
//  Version: $Id: taboos.h 856 2008-11-14 17:00:23Z davezac $
//
#ifndef _LADA_GA_TABOO_CONTAINER_H_
#define _LADA_GA_TABOO_CONTAINER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include <print/xmg.h>

namespace LaDa
{
  namespace GA
  {
    //! Holds general stuff for %Taboo operators.
    namespace Taboo
    {
      //! Creates a container which sequentially calls all of its functors.
      template< class T_INDIVIDUAL >
        class Container
        {
            //! Type of the functor container.
            typedef std::vector< boost::function<bool(const T_INDIVIDUAL&) > > t_Functors;

          public:
            //! Type of the indidividual
            typedef T_INDIVIDUAL t_Individual;
     
            //! Constructor.
            Container() {}
            //! Copy Constructor.
            Container( const Container& _c ) : functors_( _c.functors_ ) {}
            //! Virtual destrcutor, just to make sure.
            virtual ~Container() {}
     
            //! Functor over populator. Branches to correct format for t_Derived.
            bool operator()( const t_Individual& _indiv ) const
            {
              if( empty() ) return false;
              typename t_Functors :: iterator i_functor = functors_->begin();
              typename t_Functors :: iterator i_functor_end = functors_->end();
              for(; i_functor != i_functor_end; ++i_functor )
                if( (*i_functor)( _indiv ) ) return true; 
              return false;
            }
            //! Connects the callback.
            template< class T_FUNCTOR >
              void connect( const T_FUNCTOR& _functor )
              {
                if( empty() ) functors_.reset( new t_Functors );
                functors_->push_back( typename t_Functors :: value_type( _functor ) ); 
              }

            //! Returns true if container is empty.
            bool empty() const { return functors_.get() == NULL; }
            //! Clears container.
            void clear() { if( not empty() ) functors_->clear(); }

          protected:
            //! The functors.
            boost::shared_ptr<t_Functors> functors_;
        };

      namespace Factory
      {
        //! Creates a container from XML input \a _node.
        template< class T_FACTORY, class T_INDIVIDUAL >
          void create_container( T_FACTORY &_factory,
                                 const TiXmlElement &_node,
                                 Taboo::Container<T_INDIVIDUAL>& _container )
          {
            const TiXmlElement *child = _node.FirstChildElement();
            Print::xmg << Print::Xmg::comment << "Taboos" << Print::endl << Print::Xmg::indent;
            size_t nbcreated(0);
            for(; child; child = child->NextSiblingElement() )
            {
              const std::string key( child->Value() );
              if( not _factory.exists( key ) ) continue;

              boost::function<bool (const T_INDIVIDUAL&)> function;
              _factory( key, function, *child );

              if( not function.empty() ) _container.connect( function );

              ++nbcreated;
            }
            if( nbcreated == 0 )
            {
              Print::xmg << Print::Xmg::removelast << Print::Xmg::unindent;
              return;
            }
            Print::xmg << Print::Xmg :: unindent;
          }

        //! Factory for creating a container of taboos into a known taboo container.
        template< class T_FACTORY, class T_INDIVIDUAL >
          void create_container( T_FACTORY &_factory,
                                 boost::function<bool( const T_INDIVIDUAL& )>& _function,
                                 const TiXmlElement &_node,
                                 Taboo::Container<T_INDIVIDUAL>& _container )
          {
            create_container( _factory, _node, _container );
            _function = _container;
          }
        //! Factory for creating a container of taboos.
        template< class T_FACTORY, class T_INDIVIDUAL >
          void container( T_FACTORY &_factory,
                          boost::function<bool( const T_INDIVIDUAL& )>& _function,
                          const TiXmlElement &_node )
          {
            Taboo::Container<T_INDIVIDUAL> result;
            create_container( _factory, _function, _node, result );
          }
      }
    } // namespace Taboo
  } // namespace GA
} // namespace LaDa
#endif
