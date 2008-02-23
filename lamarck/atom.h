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
#include <opt/traits.h>
#include <opt/fuzzy.h>

#include <atat/vectmac.h>
#include <atat/xtalutil.h>
#include <atat/machdep.h>
#include <atat/fxvector.h>

namespace Ising_CE {

  //! \brief Describes an atom.
  //! \details An atom consists of a position, a type, and frozen status
  //!          variable. The position should always be in cartesian units. The
  //!          type can be anything, from a string with the symbol of the atom,
  //!          to an double wich codes for the atomic type somehow, to a vector
  //!          of strings which describe the possible occupations of the atomic
  //!          position. To this end, the type is a template type \a T_TYPE.
  //!          The frozen status variable indicate which, if any, coordinate
  //!          should not be touched deuring optimization. There are a number
  //!          of possibilities:
  //!            - Atom_Type::FREEZE_NONE indicate that no coordinate is
  //!                                     frozen.
  //!            - Atom_Type::FREEZE_X indicate that the cartesian x coordinate
  //!                                  is frozen.
  //!            - Atom_Type::FREEZE_Y indicate that the cartesian y coordinate
  //!                                  is frozen.
  //!            - Atom_Type::FREEZE_Z indicate that the cartesian z coordinate
  //!                                  is frozen.
  //!            - Atom_Type::FREEZE_T indicate that the occupation is frozen.
  //!            - Any combination of the above.
  //!            .
  //! \warning The default equality comparison operator compares positions only (not
  //!          occupation or site ).
  //! \xmlinput Most atoms with a simple \a T_Type can load themselves. See
  //!           \ref TagAtom.
  template<class T_TYPE>
  class Atom_Type
  {
    public:
      //! The type of the occupation
      typedef T_TYPE t_Type;

    public:
      //! No coordinat is frozen
      const static types::t_unsigned FREEZE_NONE = 0;
      //! Cartesian X coordinate is frozen ( Atom_Type::pos[0] )
      const static types::t_unsigned FREEZE_X    = 1;
      //! Cartesian Y coordinate is frozen ( Atom_Type::pos[1] )
      const static types::t_unsigned FREEZE_Y    = 2;
      //! Cartesian z coordinate is frozen ( Atom_Type::pos[2] )
      const static types::t_unsigned FREEZE_Z    = 4;
      //! Occupation coordinate is frozen ( Atom_Type::type )
      const static types::t_unsigned FREEZE_T    = 8;

    public:
      //! The atomic position in cartesian coordinate.
      atat::rVector3d pos;
      //! The atomic occupation.
      t_Type  type;
      //! The frozen status
      types::t_unsigned freeze;
      //! The site index for atoms in a structure
      types::t_int site;
      
      //! Constructor
      Atom_Type() : pos(atat::rVector3d(0,0,0)),
                    freeze(FREEZE_NONE), site(-1) {};
      //! Constructor and Initializer
      explicit 
        Atom_Type   ( const atat::rVector3d &_pos, t_Type _type) 
                  : pos(_pos), type(_type), freeze(FREEZE_NONE),
                    site(-1) {};
      //! Copy Constructor
      Atom_Type   ( const Ising_CE::Atom_Type<t_Type> &_atom )
                : pos(_atom.pos), type(_atom.type),
                  freeze(_atom.freeze), site(_atom.site) {};
      //! Equality operator. Returns true if the \e positions are equal.
      template <class TTYPE>
      bool operator== (const Atom_Type<TTYPE> &_atom) const;
      //! Compares both occupation and type.
      bool equal (const Atom_Type<t_Type> &_atom) const;
      //! Returns a constant reference to the atomic position.
      operator const atat::rVector3d& () const { return pos; }
      //! Returns a constant reference to the atomic occupation.
      operator const t_Type& () const { return type; } 
      //! Sets the atomic position.
      void operator=( const atat::rVector3d &_pos ) { pos = _pos; }
      //! Sets the occupation.
      void operator=( const t_Type &_type ) { type = _type; }

      //! Prints the atom to a stream.
      std::ostream& print_out( std::ostream &stream ) const;
      //! Loads an atom from XML.
      bool Load( const TiXmlElement &_element );
      //! Saves an atom in XML format.
      void print_xml( TiXmlElement* const node ) const;
      //! Compares the position of two atoms.
      template< class TTYPE > bool operator < ( const Atom_Type<TTYPE> &_atom ) const;
  };

  //! An atom with a complex occupation variable (eg kspace vector).
  typedef Atom_Type<std::complex<types::t_real> > CAtom;
  //! An atom with a string occupation variable.
  typedef Atom_Type<std::string> StrAtom;
  //! An atom with a double occupation variable.
  typedef Atom_Type<types::t_real> Atom;



  template<class T_TYPE> template< class TTYPE >
    inline bool Atom_Type<T_TYPE> :: operator== (const Atom_Type<TTYPE> &_atom) const
    {
      return     Fuzzy::eq( pos[0], _atom.pos[0] ) 
             and Fuzzy::eq( pos[1], _atom.pos[1] ) 
             and Fuzzy::eq( pos[2], _atom.pos[2] ); 
    }

  template<class T_TYPE> 
    inline bool Atom_Type<T_TYPE> :: equal (const Atom_Type<t_Type> &_atom) const
      {
        return     Fuzzy::eq(_atom.pos[0], pos[0] )
               and Fuzzy::eq(_atom.pos[1], pos[1] )
               and Fuzzy::eq(_atom.pos[2], pos[2] )
               and type == _atom.type;
      };

  template<class T_TYPE>
    inline std::ostream& Atom_Type<T_TYPE> :: print_out( std::ostream &stream ) const
    {
      stream << std::fixed << std::setprecision(5);
      stream << std::setw(8) << pos[0] << "  "
             << std::setw(8) << pos[1] << "  " 
             << std::setw(8) << pos[2];
      stream << "  type: " << std::setw(16)
             << Traits::Quantity< T_TYPE > :: print(type);
      if ( site != -1 )
        stream << "  site: " << site;
      stream << "  freeze: " << freeze;
      return stream;
    }

  template<class T_TYPE>
    bool Atom_Type<T_TYPE> :: Load( const TiXmlElement &_element )
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
      {
        int i=0;
        _element.Attribute("site", &i);
        site = (types::t_int) i;
      }

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
  
  template<class T_TYPE>
    void Atom_Type<T_TYPE> :: print_xml( TiXmlElement* const node ) const
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

  template<class T_TYPE> template< class TTYPE >
    bool Atom_Type<T_TYPE> :: operator < ( const Atom_Type<TTYPE> &_atom ) const
    {
      types::t_real norma = atat::norm( pos ), normb = atat::norm(_atom.pos);

      if ( norma < atat::zero_tolerance and normb > atat::zero_tolerance )
        return true;
      if ( norma < atat::zero_tolerance or normb < atat::zero_tolerance )
        return false;
       
      types::t_real a = pos[0] / norma, b = _atom.pos[0] / normb;
      if ( Fuzzy::eq(a, b) )
        return a < b;

      a = pos[1] / norma; b = _atom.pos[1] / normb;
      if ( Fuzzy::eq(a, b) )
        return a < b;

      a = pos[2] / norma; b = _atom.pos[2] / normb;
      if ( Fuzzy::eq(a, b) )
        return false;

      return a < b;
    }




//! Prints an atom to a stream.
template<class T_TYPE>
  std::ostream& operator<<(std::ostream& _stream, const Ising_CE::Atom_Type<T_TYPE> &_at )
    { _at.print_out( _stream ); return _stream; }




} // namespace Ising_CE

  
#endif
