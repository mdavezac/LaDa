//
//  Version: $Id$
//
#ifndef _ZERO_HULL_H_
#define _ZERO_HULL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <limits.h>

#include <tinyxml/tinyxml.h>

#include "opt/types.h"
#include "ch_template.h"
#include "structure.h"

namespace VA_CE {

  class One_Point_Hull : public CH_Template
  {
    public:
      const static types::t_int PRINT_XMGRACE;
      const static types::t_int PRINT_STRUCTURES;
      const static types::t_real ZERO_TOLERANCE;

    private:
      Ising_CE::Structure structure;
      types::t_real energy;

    public:
      One_Point_Hull(){ energy = (types::t_real)UINT_MAX; };
      virtual ~One_Point_Hull(){};

      // checks if structure is a breaking point.
      // if it is rebuilds the convex hull
      virtual bool add_structure( const types::t_real _energy, const Ising_CE::Structure &_structure )
      {
        if ( energy <= _energy )
          return false;
        
        // avoids redoing ch when symmetric equivalent structure is
        // found
        if ( fabs( energy - _energy ) < ZERO_TOLERANCE 
             and structure.atoms.size() <= _structure.atoms.size() )
          return false;

        energy = _energy; 
        structure = _structure;
        
        return true; 
      };


      // required lada::Fitness_Function behaviors
      virtual types::t_real evaluate() const {return 0.0;};
      virtual types::t_real evaluate(const types::t_real _x)
        { return 0.0; };
      virtual types::t_real evaluate_gradient() { return 0.0; };
      virtual types::t_real evaluate_gradient(const types::t_real _x)
        { return 0.0; };

      // other
      virtual bool Load(const TiXmlElement &_element)
      { 
        const TiXmlElement *child = _element.FirstChildElement( "BreakPoint" );
        if (not child) 
          return false;
        child->Attribute("E", &energy);
        structure.Load(*child);

        return true; 
      };
      virtual void print_out (std::ostream &stream, types::t_int print_what) const
      {
        stream << ( structure.get_concentration() )
               << " " << energy << std::endl
               << " & " << std::endl;
      }

      virtual void print_xml( TiXmlElement* const node ) const
      {
        TiXmlElement *parent;
        parent = new TiXmlElement("BreakPoint");
        parent->SetDoubleAttribute( "E", energy );
        parent->SetDoubleAttribute( "x", structure.get_concentration() );
        structure.print_xml( parent );
        node->LinkEndChild(parent);
      }
  };
  
  const types::t_int One_Point_Hull::PRINT_XMGRACE     = 1;
  const types::t_int One_Point_Hull::PRINT_STRUCTURES  = 2; 
  const types::t_real One_Point_Hull::ZERO_TOLERANCE = 1e-6; 


} // namespace VA_CE
#endif
