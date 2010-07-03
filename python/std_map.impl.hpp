#ifndef _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_
# define _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ 0
# define _PV_ ;
#endif 
#if _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 0 || _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
# ifndef _PV_ 
#   define _PV_ 
# endif

  template< class T_MAP > T_MAP* default_constructor() _PV_
#   if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
    { return new T_MAP; }
#   endif

   template< class T_MAP > T_MAP* copy_constructor( const T_MAP& _ob ) _PV_
#    if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
     { return new T_MAP(_ob); }
#    endif
     
   template< class T_MAP > T_MAP* dict_constructor( const boost::python::dict& _dict ) _PV_
#    if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
   {
     namespace bp = boost::python;
     typedef typename T_MAP :: key_type t_Key;
     typedef typename T_MAP :: value_type :: second_type t_Data;
     const size_t n( bp::len( _dict ) );
     const bp::list keys = _dict.keys();
     T_MAP* result = new T_MAP;
     for( size_t i(0); i < n; ++i )
       (*result)[ bp::extract<t_Key>(keys[i]) ] = bp::extract<t_Data>( _dict[keys[i]] );
     return result;
   }
#    endif

   template< class T_MAP > std::string print( const T_MAP& _ob ) _PV_
#    if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
   {
     std::ostringstream sstr;
     typename T_MAP::const_iterator i_var = _ob.begin();
     typename T_MAP::const_iterator i_var_end = _ob.end();
     for(; i_var != i_var_end; ++i_var )
       sstr << "(" << i_var->first << ", " << i_var->second << ") ";
     sstr << "\n";
     return sstr.str();
   }
#    endif

   template< class T_MAP > T_MAP* shallow_copy( T_MAP& _ob ) _PV_
#    if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
     { return &_ob; }
#    endif

    template< class T_MAP >
      const typename T_MAP::value_type::second_type
        getitem( const T_MAP& _ob,
                 const typename T_MAP::key_type& _k ) _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
      {
        typename T_MAP::const_iterator i_found = _ob.find( _k );
        if( _ob.end() == i_found )
        {
          static typename T_MAP::value_type::second_type none;
          std::ostringstream sstr; 
          sstr << "Key " << _k << " does not exist.\n";
          PyErr_SetString( PyExc_KeyError, sstr.str().c_str() );
          return none;
        }
        return i_found->second;
      }
#       endif


    template< class T_MAP >
      void setitem( T_MAP& _ob,
                    const typename T_MAP::key_type& _k,
                    const typename T_MAP::value_type::second_type& _d )  _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
      {  _ob[_k] = _d; }
#       endif

    template< class T_MAP >
      boost::python::object get( const T_MAP& _ob,
                                 const typename T_MAP::key_type& _k ) _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
      {
        typename T_MAP::const_iterator i_found = _ob.find( _k );
        if( _ob.end() == i_found ) return boost::python::object();
        return boost::python::object(i_found->second);
      }
#       endif

    template< class T_MAP >
      boost::python::object get_default( const T_MAP& _ob,
                                         const typename T_MAP::key_type& _k,
                                         const typename T_MAP::value_type::second_type& _d ) _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
      {
        typename T_MAP::const_iterator i_found = _ob.find( _k );
        if( _ob.end() == i_found )
          return boost::python::object(_d);
        return boost::python::object(i_found->second);
      }
#       endif

    template< class T_MAP > 
      boost::python::object set_default( T_MAP& _ob,
                                         const typename T_MAP::key_type& _k,
                                         const typename T_MAP::value_type::second_type& _d ) _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
      {
        typename T_MAP::iterator i_found = _ob.find( _k );
        if( _ob.end() == i_found )
        {
          _ob[_k] = _d;
          return boost::python::object(_ob[_k]);
        }
        return boost::python::object(i_found->second);
      }
#       endif

    template< class T_MAP > 
      boost::python::object set_default_none( T_MAP& _ob,
                                              const typename T_MAP::key_type& _k,
                                              const typename T_MAP::value_type::second_type& _d ) _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
      {
        typename T_MAP::iterator i_found = _ob.find( _k );
        if( _ob.end() == i_found )
        {
          _ob[_k] = typename T_MAP::value_type::second_type();
          return boost::python::object(_ob[_k]);
        }
        return boost::python::object(i_found->second);
      }
#       endif

    template< class T_MAP > bool has_key( T_MAP& _ob, const typename T_MAP::key_type& _k ) _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
      { return _ob.end() == _ob.find( _k ); }
#       endif


    template< class T_MAP > boost::python::list* items( T_MAP& _map ) _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
    {
      boost::python::list *list = new boost::python::list;
      typename T_MAP :: iterator i_var = _map.begin();
      typename T_MAP :: iterator i_var_end = _map.end();
      for(; i_var != i_var_end; ++i_var )
        list->append( boost::python::make_tuple( &i_var->first, &i_var->second ) );
      return list;
    }
#       endif

    template< class T_MAP > boost::python::list* values( T_MAP& _map ) _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
    {
      boost::python::list *list = new boost::python::list;
      typename T_MAP :: iterator i_var = _map.begin();
      typename T_MAP :: iterator i_var_end = _map.end();
      for(; i_var != i_var_end; ++i_var )
        list->append( &i_var->second );
      return list;
    }
#       endif

    template< class T_MAP > boost::python::list keys( T_MAP& _map ) _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
    {
      boost::python::list list; // = new boost::python::list;
      typename T_MAP :: iterator i_var = _map.begin();
      typename T_MAP :: iterator i_var_end = _map.end();
      for(; i_var != i_var_end; ++i_var )
        list.append( i_var->first.c_str() );
      return list;
    }
#       endif

    template< class T_MAP > 
      typename T_MAP :: value_type::second_type pop( T_MAP& _map,
                                                     const typename T_MAP::key_type& _k ) _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
      {
        typename T_MAP :: iterator i_var = _map.find( _k );
        if( _map.end() == i_var )
        {
          std::ostringstream sstr; 
          sstr << "Key " << _k << " does not exist.\n";
          PyErr_SetString( PyExc_KeyError, sstr.str().c_str() );
          return typename T_MAP::value_type::second_type();
        }
        typename T_MAP :: value_type::second_type  result( i_var->second );
        _map.erase( i_var );
        return result;
      }
#       endif

    template< class T_MAP > 
      typename T_MAP :: value_type::second_type popitem( T_MAP& _map ) _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
      {
        typename T_MAP :: iterator i_var = _map.begin();
        if( _map.end() == i_var )
        {
          std::ostringstream sstr; 
          sstr << "dictionary is empty.\n";
          PyErr_SetString( PyExc_KeyError, sstr.str().c_str() );
          return typename T_MAP::value_type::second_type();
        }
        typename T_MAP :: value_type::second_type result( i_var->second );
        _map.erase( i_var );
        return result;
      }
#       endif

    template< class T_MAP, template<class> class T_DEREF >
      class map_iterator _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
        : protected T_DEREF<typename T_MAP::iterator>
      {
        public:
          typedef T_MAP t_Map;
          typedef T_DEREF< typename t_Map::iterator > t_Deref;
          map_iterator() {}
          map_iterator( t_Map& _map ) : i_current( _map.begin() ), i_end( _map.end() ) {}
          map_iterator( const map_iterator& _c ) : i_current( _c.i_current ), i_end( _c.i_end ) {}
          typename t_Deref :: t_value next()
          {
            if( i_current == i_end ) 
            {
              PyErr_SetObject(PyExc_StopIteration, Py_None);
              return NULL;
            } 
            ++i_current;
            return deref( i_current );
          }
          map_iterator<T_MAP, T_DEREF> * iter() { return this; }
          bool check() const { return true; }
        protected:
          typename t_Map :: iterator i_current;
          typename t_Map :: iterator i_end;
      };
#       endif

    template< class T_ITERATOR >
      struct deref_item _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
      {
        typedef boost::python::tuple t_value;
        t_value* deref( const T_ITERATOR& _iterator )
        {
          return new boost::python::tuple( boost::python::make_tuple( _iterator->first, &_iterator->second ) );
        }
      };
#       endif
    template< class T_ITERATOR >
      struct deref_str_key _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
      {
        typedef boost::python::str t_value;
        t_value* deref( const T_ITERATOR& _iterator )
        {
          return new boost::python::str( _iterator->first.c_str() );
        }
      };
#       endif

    template< class T_ITERATOR > 
      T_ITERATOR* iter( typename T_ITERATOR :: t_Map& _map ) _PV_
#       if  _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
      { 
        return new T_ITERATOR( _map );
      }
#       endif

# ifdef _PV_ 
#   undef _PV_
# endif
# if _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 0 
#   undef _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_
#   define _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ 1
# elif _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ == 1
#   undef _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_
#   define _LADA_OPT_PYTHON_STD_MAP_IMPL_HPP_ 2
# endif
                                            
#endif
