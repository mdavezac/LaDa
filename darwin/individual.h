//
//  Version: $Id$
//
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

#include "objective.h"
#include "gatraits.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif 


namespace Individual
{
  template<class T_OBJECT, class T_QUANTITY = typename T_OBJECT :: t_Quantity,
           class T_QUANTITYTRAITS = Traits::Quantity<T_QUANTITY> >
  class Base : public eoObject, public eoPersistent
  {
    public:
      typedef T_OBJECT t_Object;
      typedef T_QUANTITYTRAITS t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_ScalarQuantity     t_ScalarQuantity;
      typedef typename t_QuantityTraits :: t_Quantity           t_Quantity;  

      typedef darwin::Fitness t_Fitness;
      typedef darwin::Fitness Fitness; // for eo


    private:
      typedef Base<t_Object, t_Quantity, t_QuantityTraits> t_Base;

    protected:
      t_Object object;
      t_Fitness repFitness;
      types::t_unsigned age;
      t_Quantity quantity;

    public: 
      Base() : age(0) {}
      
      ~Base() {}

      Base   (const t_Base &_indiv )
           : object( _indiv.object ), 
             repFitness(_indiv.repFitness),
             age(_indiv.age), quantity(_indiv.quantity) {}

      bool invalid() const
        { return repFitness.invalid(); }
      void invalidate() { repFitness.invalidate(); }

      // doesn't copy age!
      void clone( const t_Base &_indiv )
      {
        quantity   = _indiv.quantity;
        object     = _indiv.object;
        repFitness = _indiv.repFitness;
      }

      void  set_age( types::t_unsigned _age )
        { age = _age; }
      types::t_unsigned  get_age() const
        { return age; }

      bool operator<(const t_Base& _eo2) const
        { return repFitness < _eo2.repFitness; }
      bool operator>(const t_Base& _eo2) const
        { return repFitness > _eo2.repFitness; }
      bool operator==( const t_Base &_indiv ) const
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
         if ( invalid() )
           throw std::runtime_error("invalid fitness\n");
         return repFitness;
      }
      void set_fitness( types::t_real _fit )
        { repFitness = _fit; }

      // eo stuff
      virtual std::string className() const
        { return "Individual::Base<...>"; }
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


      types::t_real get_concentration() const
        { return object.get_concentration(); }

      t_Quantity& quantities() { return quantity; }
      const t_Quantity& quantities() const { return quantity; }

      t_ScalarQuantity& quantities( types::t_unsigned _n )
        { return t_QuantityTraits::scalar(quantity,_n); }
      const t_ScalarQuantity& quantities( types::t_unsigned _n ) const
        { return t_QuantityTraits::scalar(quantity,_n); }

      template<class SaveOp>
      bool Save( TiXmlElement &_node, SaveOp &_saveop ) const
      {
        TiXmlElement *xmlindiv = new TiXmlElement("Individual");
        if ( not xmlindiv ) return false;
        repFitness.Save(*xmlindiv);
        if( not _saveop(*this, *xmlindiv) ) return false;
        _node.LinkEndChild( xmlindiv );
        return true;
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

        if ( not repFitness.Load(*parent) )
          return false;

        return _loadop(*this,*parent);
      }

#ifdef _MPI
      bool broadcast( mpi::BroadCast &_bc )
      {
        return      _bc.serialize( age ) 
                and _bc.serialize( repFitness )
                and _bc.serialize<t_Object>( object )
                and t_QuantityTraits::broadcast( quantity, _bc );
      }
#endif
  };

  template<class T_OBJECT, class T_QUANTITY, class T_QUANTITYTRAITS>
  std::ostream& operator<< ( std::ostream &_str, const Base<T_OBJECT, T_QUANTITY, T_QUANTITYTRAITS> &_indiv )
  {
    typedef typename Base<T_OBJECT, T_QUANTITY, T_QUANTITYTRAITS> :: t_Object t_Object;
    std::string str; str << (const t_Object&) _indiv;
    _str << str;
    return _str;
  }

  template<class T_OBJECT>
    struct Types
    {
      typedef Base<T_OBJECT, types::t_real> Scalar; 
      typedef Base<T_OBJECT, std::vector<types::t_real> > Vector; 
    };

} // endif MultiObj

#endif
