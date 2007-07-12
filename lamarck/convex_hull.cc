#include<iomanip>

#include "convex_hull.h"
#include "atat/vectmac.h"

namespace VA_CE
{
   const types::t_int Convex_Hull::PRINT_XMGRACE     = 1;
   const types::t_int Convex_Hull::PRINT_STRUCTURES  = 2; 

  
   // adds without checking for convexity
   // if structure with same x is found, replaces
   void Convex_Hull :: force_add( const Breaking_Point &bp )
   {
     if ( break_points->empty() )
     {
       break_points->push_back(bp);
       return;
     }

     std::list<Breaking_Point> :: iterator i_begin = break_points->begin();
     std::list<Breaking_Point> :: iterator i_end = break_points->end();
     std::list<Breaking_Point> :: iterator i_insert;
     std::const_mem_fun1_ref_t<bool, Breaking_Point, const types::t_real> func( &Breaking_Point::x_geq );

     i_insert = std::find_if( i_begin, i_end, std::bind2nd(func, bp.x) );
     if ( i_insert == i_end )
       break_points->push_back(bp);
     else if ( i_insert->x == bp.x )
       *i_insert = bp;
     else 
       break_points->insert( i_insert, bp );

     while( not weed_out() );
     build_function();
   }

   
   // adds bp to convex hull, if is a new breaking point
   bool Convex_Hull :: add_structure( const Breaking_Point &bp )
   {
     std::list<Breaking_Point> :: iterator i_begin = break_points->begin();
     std::list<Breaking_Point> :: iterator i_end = break_points->end();
     std::list<Breaking_Point> :: iterator i_insert;
     std::const_mem_fun1_ref_t<bool, Breaking_Point, const types::t_real> func( &Breaking_Point::x_geq );

     // does not change convex hull if this is true
     if ( bp.E - evaluate(bp.x) > 0.0 )
       return false; 
     
     // first finds where to insert new breaking point
     i_insert = std::find_if( i_begin, i_end, std::bind2nd(func, bp.x) );
     if ( i_insert == i_end or  i_insert == i_begin ) // outside of range
       return false;
     else if( std::abs(i_insert->x - bp.x) < atat::zero_tolerance ) 
     {  // exception case, concentration already exists
       if ( std::abs(bp.E - i_insert->E) < atat::zero_tolerance )
       {
         if ( i_insert->structure.atoms.size() > bp.structure.atoms.size() )
         {
           *i_insert = bp;
           build_function();
         }
         return false;
       }

       *i_insert = bp;
     }
     else // otherwise, inserts new potential break point
       break_points->insert(i_insert, bp );


     // makes convex-hull convex
     while ( !weed_out() );
     
     // rebuilds our step-wise function
     build_function();

     return true;
   }

   // builds linear segment convex_hull
   void Convex_Hull :: build_function ()
   {
     #ifdef _DEBUG_LADA_
       if ( break_points->empty() )
       {
         std::cerr << "break_points should always include endpoints "
                   << "in void Convex_Hull :: build_function ( )"
                   << std::endl;
         exit(0);
       }
     #endif // _DEBUG_LADA_

     segments.clear(); // builds from scratch
       
     std::list<Breaking_Point> :: iterator i_bp = break_points->begin();
     std::list<Breaking_Point> :: iterator i_bp_m = i_bp;
     std::list<Breaking_Point> :: iterator i_bp_end = break_points->end();

     for ( ++i_bp; i_bp != i_bp_end; ++i_bp, ++i_bp_m )
       segments.push_back( Linear_Segment(*(i_bp_m), *i_bp) );
   }


   types::t_real Convex_Hull :: evaluate() const
   {
     #ifdef _DEBUG_LADA_
       if ( break_points->empty() )
       {
         std::cerr << " Breakpoints uninitialized in types::t_real Convex_Hull :: evaluate()" 
                   << std::endl;
         exit(0);
       }
       if ( segments.size() == 0 )
       {
         std::cerr << " Convex Hull uninitialized in types::t_real Convex_Hull :: evaluate()" 
                   << std::endl;
         exit(0);
       }
     #endif // _DEBUG_LADA_

    
     std::vector<Linear_Segment> :: const_iterator i_begin = segments.begin();
     std::vector<Linear_Segment> :: const_iterator i_end = segments.end();
     std::vector<Linear_Segment> :: const_iterator i_which;
     std::const_mem_fun1_ref_t<bool, Linear_Segment, const types::t_real> 
                              func( &Linear_Segment :: x_greater );
     i_which = std::find_if( i_begin, i_end, std::bind2nd( func, x ) );


     if ( i_which == i_end ) // exception case ...
       i_which --;

     return i_which->evaluate(x);
   }

   types::t_real Convex_Hull :: evaluate_gradient()
   {
     #ifdef _DEBUG_LADA_
       if ( break_points->empty() )
       {
         std::cerr << " Breakpoints uninitialized in " 
                   << "types::t_real Convex_Hull :: evaluate_gradient()" 
                   << std::endl;
         exit(0);
       }
       if ( segments.size() == 0 )
       {
         std::cerr << " Convex Hull uninitialized in "
                   << "types::t_real Convex_Hull :: evaluate_gradient()" 
                   << std::endl;
         exit(0);
       }
     #endif // _DEBUG_LADA_

     std::vector<Linear_Segment> :: iterator i_begin = segments.begin();
     std::vector<Linear_Segment> :: iterator i_end = segments.end();
     std::vector<Linear_Segment> :: iterator i_which;
     std::const_mem_fun1_ref_t<bool, Linear_Segment, const types::t_real> 
                              func( &Linear_Segment :: x_greater );
     i_which = std::find_if( i_begin, i_end, std::bind2nd( func, x ) );



     if ( i_which == i_end ) // exception case ...
       i_which --; 

     return i_which->a;
   }

   void Convex_Hull :: print_out (std::ostream &stream,
                                  types::t_int print_what=PRINT_XMGRACE) const
   {
     #ifdef _DEBUG_LADA_
       if ( break_points->empty() )
       {
         std::cerr << " Breakpoints uninitialized in "
                   << "types::t_real Convex_Hull :: print_out(std::ostream&)" 
                   << std::endl;
         exit(0);
       }
     #endif // _DEBUG_LADA_

     std::list<Breaking_Point> :: const_iterator i_bp   = break_points->begin();
     std::list<Breaking_Point> :: const_iterator i_last = break_points->end();

     stream << std::fixed << std::setprecision(5) << std::setw(6);
     switch (print_what)
     {
        default:
        case PRINT_XMGRACE:
          for ( ; i_bp != i_last; ++i_bp )
            i_bp->print_out(stream);
          stream << " & " << std::endl;
          break;

        case PRINT_STRUCTURES:
          for ( ; i_bp != i_last; ++i_bp )
            i_bp->structure.print_out(stream);
          break;
     }

   }

   bool Convex_Hull :: weed_out()
   {
     std::list<Breaking_Point> :: iterator i_begin, i_end, i_mid, i_last;
     types::t_real p1, p2;

     if (break_points->size() < 3 )  // two point hull always convex
       return true;

     // the following loop checks that successive triangles are convex
     i_begin = break_points->begin();
     i_last = break_points->end();
     i_mid = i_begin; ++i_mid;
     i_end = i_mid; ++i_end;
     while ( i_end != i_last ) 
     {
       p1 = ( i_begin->E - i_mid->E ) * ( i_mid->x - i_end->x );
       p2 = ( i_mid->E - i_end->E ) * ( i_begin->x - i_mid->x );
       if ( p1 >= p2 )
         break;
       ++i_begin; ++i_mid; ++i_end;
     }
     if ( i_end == i_last )  // all triangles are convex
       return true;

     // one triangle is concave -> erase midpoint
     break_points->erase(i_mid);
     
     // we don't know wether this is convex yet -> returns false
     return false; 
   }

  bool Convex_Hull :: Load( const TiXmlElement &_element )
  {
    const TiXmlElement *child;
    Breaking_Point bp;

    bp = *(break_points->begin());
    child = _element.FirstChildElement( "BreakPoint" );
    for ( ; child; child=child->NextSiblingElement( "BreakPoint" ) )
    {
      child->Attribute("E", &bp.E);
      child->Attribute("x", &bp.x);
      bp.structure.Load(*child);
      force_add(bp);
    }

    while( !weed_out() );

    build_function();

    return true;
  }


  void Convex_Hull :: print_xml( TiXmlElement &_node ) const
  {
    TiXmlElement *parent;
    
    std::list<Breaking_Point> :: const_iterator i_bp = break_points->begin();
    std::list<Breaking_Point> :: const_iterator i_end = break_points->end();
    for( ; i_bp != i_end; ++i_bp )
    {
      parent = new TiXmlElement("BreakPoint");
      parent->SetDoubleAttribute( "E", i_bp->E );
      parent->SetDoubleAttribute( "x", i_bp->x );
      if (     fabs(i_bp->x+1.0) > atat::zero_tolerance 
           and fabs(i_bp->x-1.0) > atat::zero_tolerance )
        i_bp->structure.print_xml( *parent );
      _node.LinkEndChild(parent);
    }
  }

} // namespace lada
