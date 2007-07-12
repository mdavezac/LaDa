#ifndef _DARWIN_INDIVIDUAL_H_
#define _DARWIN_INDIVIDUAL_H_

#include <exception>
#include <iostream>
#include <stdexcept>       // std::runtime_error
#include <algorithm>
#include <functional>

#include <tinyxml/tinyxml.h>

#include <eo/eoScalarFitness.h>
#include <eo/eoObject.h>      // eoObject
#include <eo/eoPersistent.h>  // eoPersistent
#include <eo/eoOpContainer.h>  


#include "opt/opt_function_base.h"
#include "opt/types.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif 


namespace darwin
{
  template<class T_OBJECT, class T_FITNESS = eoMinimizingFitness >
  class Individual : public eoObject, public eoPersistent
  {
    public:
      typedef T_OBJECT t_Object;
      typedef T_FITNESS t_Fitness;
      typedef T_FITNESS Fitness; // for eo

    protected:
      t_Object object;
      t_Fitness repFitness;
      bool quantity_is_valid;
      types::t_real quantity;
      types::t_unsigned age;

    public: 
      Individual() : quantity_is_valid (false), age(0) {}
      
      ~Individual() {}

      Individual(const Individual<t_Object, t_Fitness> &_indiv ) :
              object( _indiv.object ), 
              repFitness(_indiv.repFitness),
              quantity_is_valid (_indiv.quantity_is_valid),
              quantity(_indiv.quantity), age(_indiv.age) {}

      void invalidate() { quantity_is_valid = false; }
      bool invalid() const
        { return not quantity_is_valid; }


      void set_quantity( types::t_real _q )
      {
        quantity = _q;
        quantity_is_valid = true;
      }
      const types::t_real& get_quantity() const
        { return quantity; }

      void  set_age( types::t_unsigned _age )
        { age = _age; }
      types::t_unsigned  get_age() const
        { return age; }

      bool operator<(const Individual<t_Object, t_Fitness>& _eo2) const
        { return fitness() < _eo2.fitness(); }
      bool operator>(const Individual<t_Object, t_Fitness>& _eo2) const
        { return fitness() > _eo2.fitness(); }
      bool operator==( const Individual<t_Object, t_Fitness> &_indiv ) const
      {
        if ( std::abs(get_concentration() - _indiv.get_concentration()) > types::tolerance )
          return false;
        if ( invalid() or _indiv.invalid() )
          return object == _indiv.object; 
        if ( std::abs( quantity - _indiv.get_quantity() ) > types::tolerance )
          return false;
        return object == _indiv.object; 
      }
        
      t_Fitness fitness() const
      {
         if (not quantity_is_valid)
           throw std::runtime_error("invalid fitness\n");
         return repFitness;
      }
      void set_fitness( types::t_real _f)
        { repFitness = _f; } 

      // eo stuff
      virtual std::string className() const
        { return "Individual<...>"; }
      void printOn(std::ostream &_os) const
      {
        { // EO stuff
          if (invalid()) {
              _os << "INVALID ";
          }
          else
          {
              _os << repFitness << ' ';
          }
        } // EO stuff
      }
      void readFrom(std::istream &_is)
      {
        { // EO stuff
          std::string fitness_str;
          types::t_int pos = _is.tellg();
          _is >> fitness_str;

          if (fitness_str == "INVALID")
          {
             quantity_is_valid = false;
          }
          else
          {
             quantity_is_valid = true;
             _is.seekg(pos); // rewind
             _is >> repFitness;
          }
        } // EO stuff
      }

      void print_out( std::ostream &_stream ) const
        { object.print_out( _stream ); }

      operator t_Object& ()
        { return object; }
      operator const t_Object& () const
        { return object; }
      t_Object& Object()
        { return object; }
      const t_Object& Object() const
        { return object; }

      template<class SaveOp>
      void Save( TiXmlElement &_node, SaveOp &_saveop ) const
      {
        TiXmlElement *xmlindiv = new TiXmlElement("Individual");
        xmlindiv->SetDoubleAttribute( "fitness", fitness() );
        xmlindiv->SetDoubleAttribute( "quantity", quantity );
        _saveop(object, *xmlindiv); 
        _node.LinkEndChild( xmlindiv );
      }
      template<class LoadOp>
      bool Load( const TiXmlElement &_node, LoadOp &_loadop ) 
      {
        const TiXmlElement *parent = &_node;
        std::string name = parent->Value();
        if ( name.compare("Individual") )
          parent = _node.FirstChildElement("Individual");
        if ( not parent )
          return false;

        if ( not parent->Attribute("fitness") )
          return false;
        double d = 0;
        parent->Attribute("fitness", &d );
        repFitness = d;

        if ( not parent->Attribute("quantity") )
          return false;
        parent->Attribute("quantity", &quantity );
        quantity_is_valid = true;

        return _loadop(object,*parent);
      }

      types::t_real get_concentration() const
        { return object.get_concentration(); }

#ifdef _MPI
      bool broadcast( mpi::BroadCast &_bc )
      {
        types::t_real fitness = 0;
        if ( _bc.get_stage() == mpi::BroadCast::COPYING_TO_HERE ) fitness = repFitness;
        if (    ( not _bc.serialize( quantity_is_valid ) ) 
             or ( not _bc.serialize( quantity ) )
             or ( not _bc.serialize( age ) )
             or ( not _bc.serialize( fitness ) )
             or ( not _bc.serialize<t_Object>( object ) )   ) return false;
        if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
          repFitness = fitness;
        return true;
      }
#endif
  };


  // redefines operators < and > with respect to FITNESS
  template <class T_FITNESS, class T_OBJECT>
  bool operator<(const Individual<T_FITNESS, T_OBJECT>& _eo1,
                 const Individual<T_FITNESS, T_OBJECT>& _eo2)
    { return _eo1.operator<(_eo2); }
  template <class T_FITNESS, class T_OBJECT>
  bool operator>(const Individual<T_FITNESS, T_OBJECT>& _eo1,
                 const Individual<T_FITNESS, T_OBJECT>& _eo2)
    { return _eo2.operator<(_eo1); }

} // endif LaDa

#endif
