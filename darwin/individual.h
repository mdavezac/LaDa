#ifndef _MULTIOB_INDIVIDUAL_H_
#define _MULTIOB_INDIVIDUAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <exception>
#include <iostream>
#include <stdexcept>       // std::runtime_error
#include <algorithm>
#include <functional>

#include <tinyxml/tinyxml.h>

#include <eo/eoObject.h>      // eoObject
#include <eo/eoPersistent.h>  // eoPersistent
#include <eo/eoOpContainer.h>  


#include "opt/opt_function_base.h"
#include "opt/types.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif 


namespace Individual
{
  template<class T_OBJECT>
  class Base : public eoObject, public eoPersistent
  {
    public:
      typedef T_OBJECT t_Object;
      typedef darwin::Fitness t_Fitness;
      typedef darwin::Fitness Fitness; // for eo

    protected:
      t_Object object;
      t_Fitness repFitness;
      types::t_unsigned age;

    public: 
      Base() : age(0) {}
      
      ~Base() {}

      Base   (const Base<t_Object, t_Fitness> &_indiv )
           : object( _indiv.object ), 
             repFitness(_indiv.repFitness),
             age(_indiv.age) {}

      bool invalid() const
        { return repFitness.invalid(); }


      void  set_age( types::t_unsigned _age )
        { age = _age; }
      types::t_unsigned  get_age() const
        { return age; }

      bool operator<(const Individual<t_Object, t_Fitness>& _eo2) const
        { return repFitness < _eo2.repFitness; }
      bool operator>(const Individual<t_Object, t_Fitness>& _eo2) const
        { return repFitness > _eo2.repFitness; }
      bool operator==( const Individual<t_Object, t_Fitness> &_indiv ) const
      {
        if ( std::abs(get_concentration() - _indiv.get_concentration()) > types::tolerance )
          return false;
        if ( invalid() or _indiv.invalid() )
          return object == _indiv.object; 
        return object == _indiv.object; 
      }
        
      t_Fitness& fitness() 
        { return repFitness; }
      const t_Fitness& fitness() const
      {
         if (not quantity_is_valid)
           throw std::runtime_error("invalid fitness\n");
         return repFitness;
      }
      void set_fitness( types::t_real _fit )
        { repFitness = _fit; }

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
      void readFrom(std::istream &_is) {}

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
        repFitness->Save(*xmlindiv);
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

        if ( not repFitness->Load(*parent) )
          return false;

        return _loadop(object,*parent);
      }

      types::t_real get_concentration() const
        { return object.get_concentration(); }

#ifdef _MPI
      bool broadcast( mpi::BroadCast &_bc )
      {
        types::t_real fitness = 0;
        if ( _bc.get_stage() == mpi::BroadCast::COPYING_TO_HERE ) fitness = repFitness;
        if (    ( not _bc.serialize( age ) )
             or ( not _bc.serialize( repfitness ) )
             or ( not _bc.serialize<t_Object>( object ) )   ) return false;
        return true;
      }
#endif
  };

  template<class T_OBJECT>
  class Single : public Base<T_OBJECT>
  {
    protected:
      typedef typename Base<T_OBJECT> t_Base;
    public:
      typedef types::t_real t_Quantity;

    public:
      using t_Base :: t_Object;
      using t_Base :: t_Fitness;
      using t_Base :: Fitness;

    protected:
      t_Quantity quantity;
    
    public:
      Single () : Base<t_Object>() {}
      Single   ( const Single<t_Object> &_single )
             : Base<t_Object>( _single ), quantity( _single.quantity ) {}
    
      // doesn't copy age!
      virtual void clone( const Individual<t_Object, t_Fitness> &_indiv )
      {
        quantity   = _indiv.quantity;
        object     = _indiv.object;
        repFitness = _indiv.repFitness;
      }

      t_Quantity& quantities() { return quantities; }
      const t_Quantity& quantities() const { return quantities; }
      t_Quantity quantities( types::t_unsigned ) { return quantities; }
      const t_Quantity quantities( const types::t_unsigned ) const { return quantities; }
#ifdef _MPI
      bool broadcast( mpi::BroadCast &_bc )
      {
        if ( not _bc.serialize( quantity ) ) return false;
        return t_Base::broadcast(_bc);
      }
#endif
  };


  template<class T_OBJECT>
  class Multi : public Base<T_OBJECT>
  {
    protected:
      typedef typename Base<T_OBJECT> t_Base;
      typedef types::t_real t_Type;
    public:
      using t_Base :: t_Object;
      using t_Base :: t_Fitness;
      using t_Base :: Fitness;
      typedef std::vector<t_Type> t_Quantity;

    protected:
      t_Quantity quantity;
    
    public:
      Multi () : Base<t_Object>() {}
      Multi   ( const Single<t_Object> &_single )
            : Base<t_Object>( _single ), quantity( _single.quantity ) {}
    
      // doesn't copy age!
      virtual void clone( const Individual<t_Object, t_Fitness> &_indiv )
      {
        quantity          = _indiv.quantity;
        object            = _indiv.object;
        repFitness        = _indiv.repFitness;
      }

      t_Quantity& quantities() { return quantities; }
      const t_Quantity& quantities() const { return quantities; }
      t_Type quantities( types::t_unsigned _n ) { return quantities[_n]; }
      const t_Type quantities( const types::t_unsigned _n ) const { return quantities[_n]; }
#ifdef _MPI
      bool broadcast( mpi::BroadCast &_bc )
      {
        if ( not _bc.serialize_container( quantity ) ) return false;
        return t_Base::broadcast(_bc);
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

} // endif MultiObj

#endif
