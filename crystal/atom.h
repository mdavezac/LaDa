#ifndef _ATOM_H_
#define _ATOM_H_

#include "LaDaConfig.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <ostream>
#include <complex>
#include <boost/serialization/access.hpp>

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#include <math/eigen.h>
#include <math/traits.h>
#include <math/fuzzy.h>

#include <math/serialize.h>
#include <load_n_save/xpr/utilities.h>
#include <load_n_save/action/enum.h>

namespace LaDa
{
  namespace Crystal {

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
      friend class boost::serialization::access;
      public:
        //! The type of the occupation
        typedef T_TYPE t_Type;

        //! Tags to freeze cell coordinates.
        enum t_FreezeAtom
        {
          FREEZE_NONE =  0, //!< Freeze no atomic coordinate.
          FREEZE_X   =  1,  //!< Freeze x coordinate.
          FREEZE_Y   =  2,  //!< Freeze y coordinate.
          FREEZE_Z   =  4,  //!< Freeze z coordinate.
          FREEZE_T   =  8,  //!< Freeze type.
          FREEZE_CARTESIANS  =  7,  //!< Freeze cartesians coordinates.
          FREEZE_ALL =  15,  //!< Freeze all.
        };

      public:
        //! The atomic position in cartesian coordinate.
        math::rVector3d pos;
        //! The atomic occupation.
        t_Type  type;
        //! The frozen status
        types::t_unsigned freeze;
        //! The site index for atoms in a structure
        types::t_int site;
        
        //! Constructor
        Atom_Type() : pos(math::rVector3d(0,0,0)),
                      freeze(FREEZE_NONE), site(-1) {};
        //! Constructor and Initializer
        explicit 
          Atom_Type   ( const math::rVector3d &_pos, t_Type _type) 
                    : pos(_pos), type(_type), freeze(FREEZE_NONE),
                      site(-1) {};
        //! Copy Constructor
        Atom_Type   ( const Crystal::Atom_Type<t_Type> &_atom )
                  : pos(_atom.pos), type(_atom.type),
                    freeze(_atom.freeze), site(_atom.site) {};
        //! Equality operator. Returns true if the \e positions are equal.
        template <class TTYPE>
        bool operator== (const Atom_Type<TTYPE> &_atom) const;
        //! Compares both occupation and type.
        bool equal (const Atom_Type<t_Type> &_atom) const;
        //! Returns a constant reference to the atomic position.
        operator const math::rVector3d& () const { return pos; }
        //! Returns a constant reference to the atomic occupation.
        operator const t_Type& () const { return type; } 
        //! Sets the atomic position.
        void operator=( const math::rVector3d &_pos ) { pos = _pos; }
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
        //! To load and save to xml-like input.
        template<class T_ARCHIVE>
          bool lns_access(T_ARCHIVE const &_ar) 
          {
            namespace lns = LaDa :: load_n_save;
            std::map<std::string, LaDa::types::t_unsigned> freeze_map;
            freeze_map["none"] = FREEZE_NONE;
            freeze_map["x"] = FREEZE_X;
            freeze_map["y"] = FREEZE_Y;
            freeze_map["z"] = FREEZE_Z;
            freeze_map["t"] = FREEZE_T;
            freeze_map["cartesian"] = FREEZE_CARTESIANS;
            freeze_map["all"] = FREEZE_ALL;
            lns::xpr::Section const section = 
                         (  lns::section("Atom") 
                 << ( lns::option("pos", lns::tag=lns::required, lns::action=pos) || (
                           lns::option("x", lns::tag=lns::required, lns::action=pos[0])
                        && lns::option("y", lns::tag=lns::required, lns::action=pos[1])
                        && lns::option("z", lns::tag=lns::required, lns::action=pos[2])
                      ) 
                    )
                 << lns::option("freeze", lns::action=lns::enum_(freeze, freeze_map),
                                     lns::default_=FREEZE_NONE)
                 << lns::option("type", lns::tag=lns::required, lns::action=type) );
            return _ar & section;
          }
      private:
        //! Serializes an atom.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version);
    };

    //! An atom with a complex occupation variable (eg kspace vector).
    typedef Atom_Type<std::complex<types::t_real> > CAtom;
    //! An atom with a string occupation variable.
    typedef Atom_Type<std::string> StrAtom;
    //! An atom with a double occupation variable.
    typedef Atom_Type<types::t_real> Atom;

    template<class T_TYPE> template< class ARCHIVE >
      void Atom_Type<T_TYPE> :: serialize( ARCHIVE & _ar, const unsigned int _version)
      {
        _ar & pos;
        _ar & type;
        _ar & freeze;
        _ar & site;
      }


    template<class T_TYPE> template< class TTYPE >
      inline bool Atom_Type<T_TYPE> :: operator== (const Atom_Type<TTYPE> &_atom) const
      {
        return     math::eq( pos[0], _atom.pos[0] ) 
               and math::eq( pos[1], _atom.pos[1] ) 
               and math::eq( pos[2], _atom.pos[2] ); 
      }

    template<class T_TYPE> 
      inline bool Atom_Type<T_TYPE> :: equal (const Atom_Type<t_Type> &_atom) const
        {
          return     math::eq(_atom.pos[0], pos[0] )
                 and math::eq(_atom.pos[1], pos[1] )
                 and math::eq(_atom.pos[2], pos[2] )
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
               << math::traits::Quantity< T_TYPE > :: print(type);
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

        pos = math::rVector3d(x,y,z);

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
        types::t_real norma = pos.norm(), normb = _atom.pos.norm();

        if ( norma < types::tolerance and normb > types::tolerance )
          return true;
        if ( norma < types::tolerance or normb < types::tolerance )
          return false;
         
        types::t_real a = pos[0] / norma, b = _atom.pos[0] / normb;
        if ( math::eq(a, b) )
          return a < b;

        a = pos[1] / norma; b = _atom.pos[1] / normb;
        if ( math::eq(a, b) )
          return a < b;

        a = pos[2] / norma; b = _atom.pos[2] / normb;
        if ( math::eq(a, b) )
          return false;

        return a < b;
      }




    //! Prints an atom to a stream.
    template<class T_TYPE>
      std::ostream& operator<<(std::ostream& _stream, const Crystal::Atom_Type<T_TYPE> &_at )
        { _at.print_out( _stream ); return _stream; }

  } // namespace Crystal
} // namespace LaDa
  
#endif
