//
//  Version: $Id$
//
#ifndef _ATOM_H_
#define _ATOM_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <sstream>
#include <iomanip>
#include <ostream>
#include <complex>

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#include <atat/vectmac.h>
#include <atat/xtalutil.h>
#include <atat/machdep.h>
#include <atat/fxvector.h>

namespace Ising_CE {

  template<class T_TYPE>
  class Atom_Type
  {
    public:
      typedef T_TYPE t_Type;

    public:
      const static types::t_unsigned FREEZE_NONE;
      const static types::t_unsigned FREEZE_X;
      const static types::t_unsigned FREEZE_Y;
      const static types::t_unsigned FREEZE_Z;
      const static types::t_unsigned FREEZE_T;

    public:
      atat::rVector3d pos;
      t_Type  type;
      types::t_unsigned freeze;
      types::t_int site;
      
      Atom_Type() : pos(atat::rVector3d(0,0,0)), freeze(FREEZE_NONE), site(-1) {};
      explicit 
        Atom_Type   ( const atat::rVector3d &_pos, t_Type _type) 
                  : pos(_pos), type(_type), freeze(FREEZE_NONE), site(-1) {};
      Atom_Type   ( const Ising_CE::Atom_Type<t_Type> &_atom )
                : pos(_atom.pos), type(_atom.type), freeze(_atom.freeze), site(_atom.site) {};
      bool operator== (const Atom_Type<t_Type> &_atom) const
       { return ( pos == _atom.pos ) ? true : false; };
      bool equal (const Atom_Type<t_Type> &_atom) const
       { return pos == _atom.pos and type == _atom.type; };
      operator const atat::rVector3d& () const { return pos; }
      operator const t_Type& () const { return type; } 
      void operator=( const atat::rVector3d &_pos )
        { pos = _pos; }
      void operator=( const t_Type &_type )
        { type = _type; }

      std::ostream& print_out( std::ostream &stream ) const
      {
        stream << std::fixed << std::setprecision(5);
        stream << std::setw(8) << pos[0] << std::setw(8) << pos[1] << std::setw(8) << pos[2];
        stream << "  type: " << std::setw(16) << type;
        if ( site != -1 )
          stream << "  site: " << site;
        stream << "  freeze: " << freeze;
        return stream;
      }
      bool Load( const TiXmlElement &_element )
      {
        types::t_real x, y, z;
        if ( not _element.Attribute("x") )
         return false; 

        std::string str( _element.Attribute("x") );
        std::istringstream ss( str );
        ss >> x;

        if ( not _element.Attribute("y") )
         return false; 
        str = _element.Attribute("y");
        ss.clear(); ss.str( str );
        ss >> y;

        if ( not _element.Attribute("z") )
         return false; 
        str = _element.Attribute("z");
        ss.clear(); ss.str( str );
        ss >> z;

        pos = atat::rVector3d(x,y,z);

        ss.clear(); ss.str("");
        if ( _element.Attribute("type") )
        {
          str = _element.Attribute("type");
          ss.str( str );
        }
        ss >> type;
        

        site = -1;
        if ( _element.Attribute("site") )
          _element.Attribute("site", &site);

        freeze = FREEZE_NONE;
        if ( _element.Attribute("freeze") )
        {
          std::string str( _element.Attribute("freeze") );
          if ( str.find("all") != std::string::npos )
            freeze |= FREEZE_X | FREEZE_Y | FREEZE_Z | FREEZE_T;
          else
          {
            if ( str.find("x") != std::string::npos )
              freeze |= FREEZE_X;
            if ( str.find("y") != std::string::npos )
              freeze |= FREEZE_Y;
            if ( str.find("z") != std::string::npos )
              freeze |= FREEZE_Z;
            if ( str.find("t") != std::string::npos )
              freeze |= FREEZE_T;
          }
        }
      
        return true;
      }
      void print_xml( TiXmlElement* const node ) const
      {
        std::ostringstream ss; ss << std::fixed << std::setprecision(5);
        ss << pos[0];  node->SetAttribute( "x", ss.str().c_str() );
        ss.str(""); ss << pos[1];  node->SetAttribute( "y", ss.str().c_str() );
        ss.str(""); ss << pos[2];  node->SetAttribute( "z", ss.str().c_str() );
        ss.str(""); ss << type;    node->SetAttribute( "type", ss.str().c_str() );
        if ( freeze != FREEZE_NONE )
        {
          ss.str("");
          if ( freeze & (FREEZE_X | FREEZE_Y | FREEZE_Z | FREEZE_T) )
            ss << "all";
          else
          {
            if ( freeze & FREEZE_X )
              ss << "x";
            if ( freeze & FREEZE_Y )
              ss << "y";
            if ( freeze & FREEZE_Z )
              ss << "z";
            if ( freeze & FREEZE_T )
              ss << "t";
          }
          node->SetAttribute( "freeze", ss.str().c_str() );
        }
        if ( site > -1 )
        { 
          ss.str(""); ss << site;   node->SetAttribute( "site", ss.str().c_str() );
        }
        
      }
      bool operator < ( const Atom_Type<t_Type> &_atom ) const
      {
        types::t_real norma = atat::norm( pos ), normb = atat::norm(_atom.pos);

        if ( norma < atat::zero_tolerance and normb > atat::zero_tolerance )
          return true;
        if ( norma < atat::zero_tolerance or normb < atat::zero_tolerance )
          return false;
         
        types::t_real a = pos[0] / norma, b = _atom.pos[0] / normb;
        if ( std::abs(a-b) > atat::zero_tolerance )
          return a < b;

        a = pos[1] / norma; b = _atom.pos[1] / normb;
        if ( std::abs(a-b) > atat::zero_tolerance )
          return a < b;

        a = pos[2] / norma; b = _atom.pos[2] / normb;
        if ( std::abs(a-b) < atat::zero_tolerance )
          return false;

        return a < b;
      }
  };
  template<typename t_Type> const types::t_unsigned Atom_Type<t_Type>::FREEZE_NONE = 0;
  template<typename t_Type> const types::t_unsigned Atom_Type<t_Type>::FREEZE_X = 1;
  template<typename t_Type> const types::t_unsigned Atom_Type<t_Type>::FREEZE_Y = 2;
  template<typename t_Type> const types::t_unsigned Atom_Type<t_Type>::FREEZE_Z = 4;
  template<typename t_Type> const types::t_unsigned Atom_Type<t_Type>::FREEZE_T = 8;

  typedef Atom_Type<std::complex<types::t_real> > CAtom;
  typedef Atom_Type<std::string> StrAtom;
  typedef Atom_Type<types::t_real> Atom;
} // namespace Ising_CE

template<class T_TYPE>
  std::ostream& operator<<(std::ostream& _stream, Ising_CE::Atom_Type<T_TYPE> &_at )
    { _at.print_out( _stream ); return _stream; }
  
#endif
