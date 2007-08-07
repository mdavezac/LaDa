#ifndef _MULTIOB_OBJECTIVE_H_
#define _MULTIOB_OBJECTIVE_H_

#include<vector>
#include<math.h>
#include<stdexcept>

#include <tinyxml/tinyxml.h>

#include "opt/types.h"
#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace Objective
{
  template< class T_TYPE, class T_CONTAINER = void >
  class Base
  {
    public: 
      typedef T_TYPE t_Type;
      typedef T_CONTAINER t_Container;
    public:
      Base() {}
      virtual ~Base() {}

      virtual t_Type operator()(t_Container& ) = 0;
      virtual bool is_valid() const = 0;
      bool is_single() const { return false; }
  };
  typedef Base< types::t_real, std::vector<types::t_real> > Multi;
  template<class T_TYPE>
  class Base<T_TYPE, void>
  {
    public: 
      typedef T_TYPE t_Type;
      typedef void t_Container;
    public:
      Base() {}
      virtual ~Base() {}


      virtual t_Type operator()(t_Type ) = 0;
      virtual bool is_valid() const = 0;
      bool is_single() const { return true; }
  };
  typedef Base<types::t_real> Single;

  class Maximize : public Single
  {
    public:
      Maximize() {}
      virtual ~Maximize() {}
      
      virtual t_Type operator()(t_Type _val) { return -_val; }
      bool is_valid() const { return true; }
  };
  class Minimize : public Single
  {
    public:
      Minimize() {}
      virtual ~Minimize() {}
      
      virtual t_Type operator()(t_Type _val) { return _val; }
      bool is_valid() const { return true; }
  };
  class Target : public Single
  {
    protected:
      t_Type target; 
    public:
      Target( t_Type _target ) : target( _target ) {}
      virtual ~Target() {}
      
      virtual t_Type operator()(t_Type _val) 
      { 
        return std::abs( _val - target );
      }
      bool is_valid() const { return true; }
  };


  class Container : public Multi
  {
    protected:
      typedef std::vector< Single* > t_objectives;

    protected:
      t_objectives objectives;

    public:
      Container() {}
      virtual ~Container () { DestroyContainer(); }

      bool is_valid() const
      {
        t_objectives::const_iterator i_objective = objectives.begin();
        t_objectives::const_iterator i_end = objectives.begin();
        for(; i_objective != i_end; ++i_objective )
          if ( not (*i_objective)->is_valid() )
            return false;
        return true;
      }
      
    private:
      void DestroyContainer()
      {
        t_objectives :: iterator i_objective = objectives.begin();
        t_objectives :: iterator i_end = objectives.end();
        for(; i_objective != i_end; ++i_objective )
          delete *i_objective;
        objectives.clear();
      }
  };

  class LinearSum : public Container
  {
    protected:
      std::vector<types::t_real> coefs;

    public:
      LinearSum() {}
      virtual ~LinearSum() {}

      virtual types::t_real operator()(t_Container& _val) 
      {
        if ( _val.size() != coefs.size() )
          throw std::runtime_error( "Wrong number of objective functions\n" );

        types::t_real inter = 0;
        t_Container::iterator i_val = _val.begin();
        std::vector< types::t_real >::iterator i_coef = coefs.begin();
        t_objectives::iterator i_objective = objectives.begin();
        t_objectives::iterator i_end = objectives.begin();
        for(; i_objective != i_end; ++i_objective, ++i_coef, ++i_val )
          inter += ( *i_coef ) * (*i_objective)->operator()( *i_val );
          
        return inter;
      };
  };

  Single *new_from_xml( const TiXmlElement &_node );
}  // namespace fitness

namespace darwin
{
  class Fitness
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<darwin::Fitness>( darwin::Fitness & ); 
#endif 
    public:
      typedef types::t_real t_Quantity;

    protected:
      t_Quantity quantity;
      bool is_valid;

    public:
      Fitness() : is_valid( false )  {}
      Fitness( const Fitness & _c ) : quantity( _c.quantity ), is_valid( false ) {}
      ~Fitness() {}


      bool operator>(const Fitness & _f)
        { return quantity > _f.quantity; }
      bool operator<(const Fitness & _f)
        { return quantity < _f.quantity; }
      bool operator==(const Fitness & _f)
        { return std::abs(quantity - _f.quantity) < types::tolerance; }

      bool invalid() const { return not is_valid; }
      void invalidate() { is_valid = false; }
      void operator=( types::t_real _fit ) 
      {
        is_valid = true;
        quantity = _fit;
      }
      operator t_Quantity() const
      { 
        if ( not is_valid )
          throw std::runtime_error( " Invalid Fitness !!\n" );
        return quantity;
      }
      bool Load( const TiXmlElement & _node );
      bool Save( TiXmlElement & _node ) const;
  };
}

std::ostream & operator<<( std::ostream &_os, darwin::Fitness &_fit );

#endif //  _MULTIOB_OBJECTIVE_H_
