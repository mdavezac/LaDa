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
#include "fitness.h"
#include "gatraits.h"


#ifdef _MPI
#include "mpi/mpi_object.h"
#endif 


namespace Individual
{
  template<class T_INDIVTRAITS>
  class Base : public eoObject, public eoPersistent
  {
    public:
      typedef T_INDIVTRAITS t_IndivTraits;
      typedef typename t_IndivTraits::t_Fitness t_Fitness;
      typedef typename t_IndivTraits::t_ScalarFitness Fitness;
    protected:
      typedef typename t_IndivTraits :: t_Object t_Object;
      typedef typename t_IndivTraits :: t_QuantityTraits t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_ScalarQuantity     t_ScalarQuantity;
      typedef typename t_QuantityTraits :: t_Quantity           t_Quantity;  

    private:
      typedef Base<t_IndivTraits> t_Base;

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
      void set_fitness( const t_Quantity _q )
        { repFitness = _q; }

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
        { _stream << (const t_Object&) object << "  "; }

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
        if ( not xmlindiv )
        {
          std::cerr << "Memory Allocation Error while saving individual" << std::endl;
          return false;
        }

        if (     repFitness.Save(*xmlindiv)
             and _saveop(*this, *xmlindiv)  )
        {
          _node.LinkEndChild( xmlindiv );
          return true;
        }

        delete xmlindiv;
        std::cerr << "Errow while save individual" << std::endl <<  *this << std::endl;
        return false;
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
                and repFitness.broadcast(_bc)
                and _bc.serialize<t_Object>( object )
                and t_QuantityTraits::broadcast( quantity, _bc );
      }
#endif
  };

  template< class T_BASE >
  class Multi : public T_BASE
  {
      typedef T_BASE t_Base;
    public: 
      typedef typename t_Base::t_IndivTraits t_IndivTraits; 

    protected:
      typedef typename t_Base::t_Fitness::t_ScalarFitness t_ScalarFitness;
      typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;

    protected:
      using t_Base :: repFitness;

    public:
      void set_fitness( const t_ScalarFitnessQuantity _q )
        { (t_ScalarFitness&) repFitness = _q; }
      template<class SaveOp>
      bool Save( TiXmlElement &_node, SaveOp &_saveop ) const
      {
        TiXmlElement *xmlindiv = new TiXmlElement("Individual");
        if ( not xmlindiv )
        {
          std::cerr << "Memory Allocation Error while saving individual" << std::endl;
          return false;
        }

        if (     repFitness.Save(*xmlindiv)
             and _saveop(*this, *xmlindiv)  )
        {
          _node.LinkEndChild( xmlindiv );
          return true;
        }

        delete xmlindiv;
        std::cerr << "Errow while save individual" << std::endl <<  *this << std::endl;
        return false;
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
  };

  template<class T_INDIVTRAITS>
  inline std::ostream& operator<< ( std::ostream &_str,
                                    const Base<T_INDIVTRAITS> &_indiv )
    { _indiv.print_out(_str); return _str; }

  template<class T_OBJECT>
    struct Types
    {
        typedef T_OBJECT t_Object;
        typedef Base< Traits :: Indiv< t_Object, Traits::Quantity<types::t_real> > > Scalar;
      protected:
        typedef Base< Traits :: Indiv< t_Object,
                                       Traits::Quantity< std::vector<types::t_real> > > > VectorBase;
      public:
        typedef Multi< VectorBase > Vector;
    };

} // endif MultiObj

#endif
