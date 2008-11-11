//
//  Version: $Id$
//
#ifndef _LADA_OPT_FACTORY_H_
#define _LADA_OPT_FACTORY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_map.hpp>

namespace LaDa
{
  //! Holds factory related objects.
  namespace Factory
  {
    //! A policy for chaining connect calls.
    template< class T_THIS > 
    class ChainConnects
    {
        friend class T_THIS;
      public:
        //! Type of the derived class.
        typedef T_THIS t_This;
       
        //! Functor which chains calls to connect.
        template< class T_FUNCTOR >
          ChainConnects& operator()( typename t_This::t_Key& _key,
                                     const T_FUNCTOR& _functor )
            { this_.connect( _key, _functor ); return *this; }

        //! Copy constructor. 
        ChainConnects( const ChainConnects &_c ) : this_( _c.this_ ) {}

      private:
        //! Private constructor. 
        ChainConnects( t_This &_this ) : this_( _this ) {}
        //! Holds a reference to PureCalls.
        t_This &this_;
    };

    
    //! \brief Makes pure calls to functors taking no argument and returning no
    //!        values. Functors must be copy constructible.
    class PureCalls : public ChainConnects< PureCalls >
    {
      //! Pure abstract class for storing purposes.
      class BaseType;
      //! The derived objects which actually store the functors.
      template< class T_FUNCTOR > class DerivedType;
      //! The type of the key.
      typedef const std::string t_Key;
      //! The type of the map.
      typedef boost::ptr_map< t_Key, BaseType > t_Map;
      
      public:
        //! Constructor.
        PureCalls() {}
        //! virtual Destructor.
        virtual ~PureCalls() {}

        //! \brief Adds a new functor. 
        //! \details Throws on duplicates.
        template< class T_FUNCTOR >
          ChainConnects<PureCalls> connect( const t_Key& _key, const T_FUNCTOR& _functor );

        //! performs the call.
        void operator()( const t_Key& _key );

        //! \brief Deletes a connection.
        //! \details Unlike other member functions, this one does not throw if
        //!          \a _key does not exist..
        void disconnect( const t_Key& _key );
         
      protected:
        //! The map.
        t_Map map_;
    };

    struct PureCalls :: BaseType 
    {
      //! The virtual functor. 
      virtual void operator()() = 0;
      //! A virtual Destructor.
      virtual ~BaseType() {};
    };

    template< class T_FUNCTOR >
      class PureCalls :: DerivedType : public BaseType 
      {
        public:
          //! the type of the functor.
          typedef T_FUNCTOR t_Functor;
          //! Constructor
          DerivedType   ( const t_Functor& _functor )
                      : functor_( new t_Functor( _functor ) ) {}
          //! Copy Constructor.
          DerivedType( const DerivedType& _c ) : functor_( _c.functor_ ) {}
          //! Destructor.
          virtual ~DerivedType() {}
         
          //! The virtual functor. 
          void operator()() { (*functor_)(); }

        protected:
          //! Holds the functor.
          boost::shared_ptr< t_Functor > functor_;
      };

    template< class T_FUNCTOR >
      ChainConnects<PureCalls> PureCalls :: connect( const t_Key& _key,
                                                     const T_FUNCTOR& _functor )
      {
        __DOASSERT( map_.end() != map_.find( _key ),
                    "Key " << _key << " already exists.\n" )
        typedef typename t_Map :: value_type value_type;
        map_.insert( _key, new DerivedType<T_FUNCTOR>( _functor ) );
        return ChainConnects<PureCalls>( *this );
      }

    inline void PureCalls :: operator()( const t_Key& _key )
    {
       t_Map :: iterator i_functor = map_.find( _key );
       __DOASSERT( i_functor == map_.end() , "Key " << _key << " does not exists.\n" )
       (*i_functor->second)();
    }

    inline void PureCalls :: disconnect( const t_Key& _key )
    {
       t_Map :: iterator i_functor = map_.find( _key );
       if( i_functor != map_.end() ) map_.erase( i_functor );
    }

  }
}

#endif 
