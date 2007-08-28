//
//  Version: $Id$
//
#ifndef _CONVEX_HULL_H_
#define _CONVEX_HULL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <list>
#include <iostream>

#include <tinyxml/tinyxml.h>

#include "opt/types.h"
#include "structure.h"
#include "ch_template.h"

namespace VA_CE {

  struct Breaking_Point
  {
    types::t_real E,x;
    Ising_CE::Structure structure;
    Breaking_Point() : E(-666.666), x(-666.666), structure() {};
    Breaking_Point   (const types::t_real _E, const Ising_CE::Structure &_structure) 
                   : E(_E), structure(_structure)
    { x = structure.get_concentration(); }
    Breaking_Point(const Breaking_Point &_bp ) : E(_bp.E), x(_bp.x), structure(_bp.structure) {};

    bool x_less(const Breaking_Point &bp) const
      { return (x < bp.x); }
    bool x_less(const types::t_real _x) const
      { return (x < _x); }
    bool x_greater(const types::t_real _x) const
      { return (x > _x); }
    bool x_geq(const types::t_real _x) const // greater or equal
      { return ( std::abs(x-_x) < atat::zero_tolerance ) ? true : (x > _x); }
    bool x_less(const Breaking_Point &bp1, const Breaking_Point &bp2) const
      { return (bp1.x < bp2.x); }
    bool E_less(const Breaking_Point &bp) const
      { return (E < bp.E); }
    bool E_less(const Breaking_Point &bp1, const Breaking_Point &bp2) const
      { return (bp1.E < bp2.E); }
    bool operator==(const Breaking_Point &bp) const
    {
      if (    std::abs(E - bp.E) > atat::zero_tolerance
           or std::abs( x - bp.x ) > atat::zero_tolerance )
        return false;
      return (structure == bp.structure); 
    }

    void operator=(const Breaking_Point &bp)
      { E = bp.E; x = bp.x; structure = bp.structure; }

    void print_out(std::ostream &stream) const
      { stream << x << " " << E << std::endl; }

    types::t_real get_concentration() const
      { return structure.get_concentration(); }
  };

  class Convex_Hull : public CH_Template
  {
    public:
      const static types::t_int PRINT_XMGRACE;
      const static types::t_int PRINT_STRUCTURES;

    protected:
      
      struct Linear_Segment
      { 
        types::t_real a, b;
        types::t_real x_end;
        types::t_real evaluate(const types::t_real x) const { return (a*x + b); }
        types::t_real evaluate_gradient() const { return a; }
        Linear_Segment() : a(0), b(0), x_end(0) {};
        Linear_Segment( const Breaking_Point &start, const Breaking_Point &end)
        {
          b = start.x - end.x; 
          if ( b == 0 ) 
            { std::cerr << "Error initialising linear segment" << std::endl; exit(0); }
          a = (start.E - end.E)/b; b = end.E - end.x*a;
          x_end = end.x;
        }
        bool x_greater(const types::t_real _x) const
          { return (_x < x_end); }
      };

    protected:
      types::t_real x; // concentration at which evaluate and evaluate_gradient are given
      std::list<Breaking_Point> *break_points;
      std::vector<Linear_Segment> segments;

    protected:
      void build_function();
      bool weed_out(); // removes unconvex points one by one

    public:
      Convex_Hull() : segments()
      {
        segments.reserve(20); 
        break_points = new std::list<Breaking_Point>;
      }
      virtual ~Convex_Hull()
      {
        if ( break_points ) 
          delete break_points;
        break_points = NULL;
      };

      // checks if structure is a breaking point.
      // if it is rebuilds the convex hull
      bool add_structure( const Breaking_Point &bp );
      virtual bool add_structure( const types::t_real E, const Ising_CE::Structure &str )
      { Breaking_Point bp(E,str); return add_structure( bp ); }
        // adds structure without checking for convexity -- dangerous
      void force_add( const Breaking_Point &bp );


      // required lada::Fitness_Function behaviors
      virtual types::t_real evaluate() const;
      virtual types::t_real evaluate(const types::t_real _x)
        { x = _x; return evaluate(); }
      virtual types::t_real evaluate_gradient();
      virtual types::t_real evaluate_gradient(const types::t_real _x)
        { x = _x; return evaluate_gradient(); }

      // other
      virtual void print_out (std::ostream &stream, types::t_int print_what) const;

      virtual bool Load(const TiXmlElement &_element);
      virtual void print_xml( TiXmlElement &_node ) const;
      void clear() { break_points->clear(); }
  };



} // namespace lada
#endif
