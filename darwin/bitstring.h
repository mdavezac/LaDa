//
//  Version: $Id$
//
#ifndef _BITSTRING_H_
#define _BITSTRING_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<string>
#include<sstream>
#include <functional>

#include<tinyxml/tinyxml.h>

#include<opt/types.h>

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace BitString
{

  template< class T_CONTAINER = std::vector< types::t_real > >
  struct Object
  {
    public:
      typedef T_CONTAINER t_Container;
      typedef typename t_Container :: value_type t_Type;
      typedef typename t_Container :: iterator iterator;
      typedef typename t_Container :: const_iterator const_iterator;

    public:
      t_Container bitstring;

    public:
      Object() {}
      Object(const Object &_c) : bitstring(_c.bitstring) {};
      Object(const t_Container &_c) : bitstring(_c) {};
      ~Object() {};
    
      bool operator==( const Object &_c ) const
      {
        return std::equal( bitstring.begin(), bitstring.end(), 
                            _c.bitstring.begin() ); 
      }

      void mask( types::t_unsigned _start, types::t_unsigned _end);
      types::t_unsigned size() { return bitstring.size(); }
      t_Container& Container() { return bitstring; }
      const t_Container& Container() const { return bitstring; }

      const_iterator begin() const
        { return bitstring.begin();  }
      const_iterator end() const
        { return bitstring.end();  }
      iterator begin() 
        { return bitstring.begin();  }
      iterator end()
        { return bitstring.end();  }
      
      types::t_real get_concentration() const
      {
        typename t_Container :: const_iterator i_var = bitstring.begin();
        typename t_Container :: const_iterator i_var_end = bitstring.end();
        types::t_real result = 0.0;
        for(; i_var != i_var_end; ++i_var )
          result += *i_var > 0 ? 1.0: -1.0;
        result /= static_cast<types::t_real>(bitstring.size());
        return result;
      }

#ifdef _MPI
       bool broadcast ( mpi::BroadCast &_bc )
       {
         return _bc.serialize_container( bitstring );
       }
#endif
  };

  template<class T_CONT>
    std::ostream& operator<<(std::ostream &_stream, const Object<T_CONT> &_o);
  template<class T_CONT>
    void operator<<(std::string &_str, const Object<T_CONT> &_o);
  template<class T_CONT>
    void operator<<(Object<T_CONT> &_o, const std::string &_str );

  template< class T_TYPE > inline bool spin_up( T_TYPE _t ) { return _t > T_TYPE(0); }
  template<> inline bool spin_up<bool>(bool _t ) { return _t; }
  template< class T_TYPE > inline bool spin_down( T_TYPE _t ) 
    { return  not spin_up<T_TYPE>(_t); }
  template<> inline bool spin_down<bool>( bool _t ) { return not _t; }

  template< class T_TYPE >
    void inline flip( T_TYPE &_t ) { _t = -_t; }
  template<>
    void inline flip< bool >( bool &_t ) { _t = not _t; }

  template< class T_CONTAINER >
    inline void Object<T_CONTAINER> :: mask( types::t_unsigned _start,
                                             types::t_unsigned _end)
    {
      if ( _end > bitstring.size() ) _end = bitstring.size();
      std::for_each( bitstring.begin() + _start, bitstring.begin() + _end,
                     std::ptr_fun( &flip< typename T_CONTAINER::value_type > ) );  
    }
    
  template< class T_OBJECT >
    class Crossover 
    {
      public:
        typedef T_OBJECT t_Object;
      public:
        types::t_real rate;

      public:
        Crossover() : rate (0.5) {}
        Crossover( types::t_real &_rate ) : rate ( _rate ) {}
        Crossover( const TiXmlElement &_node ) : rate ( 0.5 ) { Load( _node ); }

        bool Load( const TiXmlElement &_node )
        {
          _node.Attribute("rate", &rate);
          if ( rate <= types::tolerance ) rate = 0.5;
          if ( 1.0 - rate <= types::tolerance ) rate = 1.0;
          return true;
        }
        Crossover( const Crossover<t_Object> &_k ) : rate (_k.rate) {}
        ~Crossover() {}

        bool operator()( t_Object &_off, const t_Object &_parent );
        bool operator()( t_Object &_off1, t_Object &_off2 );
        std::string print_out() const
        {
          std::ostringstream sstr;
          sstr << "Crossover rate = " << rate;
          return sstr.str();
        }
        std::string className() const { return "BitString::Crossover"; }
    };

  template< class T_OBJECT >
    inline bool Crossover<T_OBJECT> :: operator()( t_Object &_off, const t_Object &_parent )
    {
      if( rate <= types::tolerance  ) return false;
      if( 1.0 - rate <= types::tolerance )
      {
        _off = _parent;
        return true;
      }

      typename t_Object :: iterator  i_var_off = _off.bitstring.begin();
      typename t_Object :: iterator  i_var_off_end = _off.bitstring.end();
      typename t_Object :: const_iterator  i_var_par = _parent.bitstring.begin();
      typename t_Object :: const_iterator  i_var_par_end = _parent.bitstring.end();
      bool has_changed = false;
      for(; i_var_off != i_var_off_end and i_var_par != i_var_par_end; ++i_var_off, ++i_var_par )
      {
        if( not rng.flip( rate ) ) continue;
        *i_var_off = *i_var_par;
        has_changed = true;
      }
      return has_changed;
    }
  template< class T_OBJECT >
    inline bool Crossover<T_OBJECT> :: operator()( t_Object &_off1, t_Object &_off2 )
    {
      if( rate <= types::tolerance )
      {
        _off2 = _off1;
        return true;
      }
      if( 1.0 - rate <= types::tolerance )
      {
        _off1 = _off2;
        return true;
      }

      typename t_Object :: iterator  i_var_1 = _off1.bitstring.begin();
      typename t_Object :: iterator  i_var_1_end = _off1.bitstring.end();
      typename t_Object :: iterator  i_var_2 = _off2.bitstring.begin();
      typename t_Object :: iterator  i_var_2_end = _off2.bitstring.end();
      for(; i_var_1 != i_var_1_end and i_var_2 != i_var_2_end; ++i_var_1, ++i_var_2 )
        if( rng.flip( rate ) ) *i_var_1 = *i_var_2;
        else                   *i_var_2 = *i_var_1;
      return true;
    }
   
  template< class T_OBJECT >
    class Mutation 
    {
      public:
        typedef T_OBJECT t_Object;
      public:
        types::t_real rate;

      public:
        Mutation() : rate (0.1) {}
        Mutation( types::t_real &_rate ) : rate ( _rate ) {}
        Mutation( const TiXmlElement &_node ) : rate ( 0.5 ) { Load(_node); }
          
        bool Load(const TiXmlElement &_node )
        {
          _node.Attribute("rate", &rate);
          if ( rate <= types::tolerance ) rate = 0.5;
          if ( 1.0 - rate <= types::tolerance ) rate = 1.0;
          return true;
        }
        Mutation( const Mutation<t_Object> &_k ) : rate (_k.rate) {}
        ~Mutation() {}

        bool operator()( t_Object &_o);
        std::string print_out() const
        {
          std::ostringstream sstr;
          sstr << "Mutation rate = " << rate;
          return sstr.str();
        }
        std::string className() const { return "BitString::Mutation"; }
    };


  template< class T_OBJECT >
    inline bool Mutation<T_OBJECT> :: operator()( t_Object &_o )
    {
      if( rate <= types::tolerance  ) return false;
      if( 1.0 - rate <= types::tolerance ) rate = 1.0;

      typename t_Object :: iterator  i_var = _o.bitstring.begin();
      typename t_Object :: iterator  i_var_end = _o.bitstring.end();
      bool has_changed = false;
      for(; i_var != i_var_end; ++i_var )
      {
        if( not rng.flip( rate ) ) continue;
        flip< typename t_Object :: t_Container :: value_type >( *i_var );
        has_changed = true;
      }
      return has_changed;
    }

}

#include "bitstring.impl.h"

#endif
