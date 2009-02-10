//
//  Version: $Id$
//
#ifndef _LADA_C2FHANDLES_H_
#define _LADA_C2FHANDLES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector> 
#include <list> 

#include <crystal/atom.h>
#include <crystal/structure.h>
#include <opt/types.h>
#include <opt/debug.h>


namespace LaDa
{
  //! Creates a handle manager to go from c to fortran.
  template< class T_OBJECT >
    class C2FHandles
    {
      public:
        //! Type of the object to handle.
        typedef T_OBJECT t_Object;

        //! Constructor and Initializer
        C2FHandles() {}
        //! Copy Constructor
        C2FHandles( const C2FHandles &_c ) {}
        //! \brief Destructor
        ~C2FHandles() {}

        //! Adds a handle/object.
        size_t push_back( const T_OBJECT &_object )
        { 
          objects().push_back( _object );
          handles().push_back( --(objects().end()) );
          const t_Objects &cobjects( objects() );
          chandles().push_back( --(cobjects.end()) );
          exists().push_back( true );
          return handles().size() - 1;
        }
        //! Erases a handle/object.
        void erase( size_t _index )
        { 
          __DOASSERT( exists().size() <= _index, "Index out of range.\n" )
          __DOASSERT( not exists()[_index], "Handle is invalid.\n" )
          objects().erase( handles()[_index] );
          exists()[_index] = false;
        }
        //! Returns an object.
        t_Object& operator[]( size_t _index )
        {
          __DOASSERT( exists().size() <= _index, "Index out of range.\n" )
          __DOASSERT( not exists()[_index], "Handle is invalid.\n" )
          return *( handles()[_index] );
        }
        //! Returns an object.
        const t_Object& operator[]( size_t _index ) const
        {
          __DOASSERT( exists().size() <= _index, "Index out of range.\n" )
          __DOASSERT( not exists()[_index], "Handle is invalid.\n" )
          return *( chandles()[_index] );
        }
      protected:
        //! Type of the container of objects.
        typedef std::list< t_Object > t_Objects;
        //! Type of the container of handles.
        typedef std::vector< typename t_Objects :: iterator > t_Handles;
        //! Type of the container of constant handles.
        typedef std::vector< typename t_Objects :: const_iterator > t_cHandles;
        //! Type of the container of existence.
        typedef std::vector< bool > t_Exists;
        //! A static singleton for the handles.
        static t_Handles& handles();
        //! A static singleton for the constant handles.
        static t_cHandles& chandles();
        //! A static singleton for the objects.
        static t_Objects& objects();
        //! A static singleton for existence.
        static t_Exists& exists();
    };

  template< class T_OBJECT > 
    typename C2FHandles<T_OBJECT> :: t_Handles& C2FHandles<T_OBJECT> :: handles()
    {
      static t_Handles container_;
      return container_;
    }
  template< class T_OBJECT > 
    typename C2FHandles<T_OBJECT> :: t_cHandles& C2FHandles<T_OBJECT> :: chandles()
    {
      static t_cHandles container_;
      return container_;
    }
  template< class T_OBJECT > 
    typename C2FHandles<T_OBJECT> :: t_Objects& C2FHandles<T_OBJECT> :: objects()
    {
      static t_Objects container_;
      return container_;
    }
  template< class T_OBJECT > 
    typename C2FHandles<T_OBJECT> :: t_Exists& C2FHandles<T_OBJECT> :: exists()
    {
      static t_Exists container_;
      return container_;
    }
} // namespace LaDa

#endif // _VFF_FUNCTIONAL_H_
