//
//  Version: $Id$
//

#include <boost/assign/ptr_map_inserter.hpp>

namespace LaDa
{
  namespace GA
  {
    namespace Factory
    {

      template< class T_POPULATOR, class T_ARG >
        struct Operators<T_POPULATOR, T_ARG> :: BaseType 
        {
          //! Type of the argument.
          typedef typename Operators<T_POPULATOR, T_ARG>::t_Function t_Function;
          //! The virtual functor. 
          virtual void create_basetype( t_Function &_function, T_ARG& _arg ) = 0;
          //! A virtual Destructor.
          virtual ~BaseType() {};
        };
      template< class T_POPULATOR, class T_ARG > template< class T_FUNCTOR >
        class Operators<T_POPULATOR, T_ARG> :: DerivedType 
            : public Operators<T_POPULATOR, T_ARG>::BaseType 
        {
            //! Type of the argument.
            typedef typename Operators<T_POPULATOR, T_ARG>::t_Function t_Function;
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
            void create_basetype( t_Function& _result, T_ARG&_arg )
              { (*functor_)( _result, _arg ); }
 
          protected:
            //! Holds the functor.
            boost::shared_ptr< t_Functor > functor_;
        };
 
      template< class T_POPULATOR, class T_ARG > template< class T_FUNCTOR >
        ::LaDa::Factory::ChainConnects< Operators<T_POPULATOR, T_ARG> >
          Operators<T_POPULATOR, T_ARG> :: connect( const t_Key& _key,
                                                    const T_FUNCTOR& _functor )
          {
            __DOASSERT( exists( _key ),
                        "Key " << _key << " already exists.\n" )
            boost::assign::ptr_map_insert< DerivedType<T_FUNCTOR> >
              ( map_ )( _key, _functor );
            return ::LaDa::Factory::ChainConnects< t_This >( *this );
          }

      template< class T_POPULATOR, class T_ARG >
        void Operators<T_POPULATOR, T_ARG> :: operator()( const t_Key& _key,
                                                          t_Function &_result,
                                                          t_Arg &_arg )
        {
           __DOASSERT( not exists( _key ) ,
                       "Key " << _key << " does not exists.\n" )
           typename t_Map :: iterator i_functor = map_.find( _key );
           (*i_functor->second).create_basetype( _result, _arg );
        }
 
      template< class T_POPULATOR, class T_ARG >
        void Operators<T_POPULATOR, T_ARG> :: disconnect( const t_Key& _key )
        {
           typename t_Map :: iterator i_functor = map_.find( _key );
           if( i_functor != map_.end() ) map_.erase( i_functor );
        }

      template< class T_POPULATOR, class T_ARG >
        std::ostream& operator<<( std::ostream& _stream, 
                                  const Operators<T_POPULATOR, T_ARG>& _factory )
        {
          typedef Operators<T_POPULATOR, T_ARG> t_Factory;
          typename t_Factory::t_Map::const_iterator i_map = _factory.map_.begin();
          typename t_Factory::t_Map::const_iterator i_map_end = _factory.map_.end();
          _stream << "Factory keys: \n";
          for(; i_map != i_map_end; ++i_map )
            _stream << "  _ " << i_map->first << "\n";
          return _stream;
        }
    }
  }
}

