#ifndef _ATOM_H_
#define _ATOM_H_

#include "LaDaConfig.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <ostream>
#include <complex>

#include <boost/serialization/access.hpp>

#include <opt/types.h>
#include <math/eigen.h>
#include <math/fuzzy.h>

#include <math/serialize.h>
# ifdef LADA_WITH_LNS
#  include <load_n_save/lns.h>
# endif

namespace LaDa
{
  namespace Crystal
  {

    //! \brief Describes an atom.
    //! \details An atom consists of a position and a type. The position should
    //!          always be in cartesian units. The type can be anything, from a
    //!          string with the symbol of the atom, to an double wich codes
    //!          for the atomic type somehow, to a vector of strings which
    //!          describe the possible occupations of the atomic position. To
    //!          this end, the type is a template type \a T_TYPE. 
    //! \warning The default equality comparison operator compares positions only (not
    //!          occupation or site ).
    template<class T_TYPE>
    class AtomBase
    {
      friend class boost::serialization::access;
#     ifdef LADA_WITH_LNS
#       friend class load_n_save::access;
#     endif
      public:
        //! The type of the occupation
        typedef T_TYPE t_Type;

      public:
        //! The atomic position in cartesian coordinate.
        math::rVector3d pos;
        //! The atomic occupation.
        t_Type  type;
        
        //! Constructor
        Atom_Type() : pos(math::rVector3d(0,0,0)), type() {};
        //! Constructor and Initializer
        explicit AtomBase   ( const math::rVector3d &_pos, t_Type _type) 
                          : pos(_pos), type(_type) {};
        //! Copy Constructor
        AtomBase(const AtomBase &_c) : pos(_c.pos), type(_c.type) {};
        //! Equality operator. Returns true if the \e positions are equal.
        template <class TTYPE> bool operator== (const AtomBase<TTYPE> &_atom) const;
        //! Compares both occupation and type.
        bool equal (const AtomBase<t_Type> &_atom) const;
        //! Returns a constant reference to the atomic position.
        operator const math::rVector3d& () const { return pos; }
        //! Returns a constant reference to the atomic occupation.
        operator const t_Type& () const { return type; } 
        //! Sets the atomic position.
        void operator=( const math::rVector3d &_pos ) { pos = _pos; }
        //! Sets the occupation.
        void operator=( const t_Type &_type ) { type = _type; }

        //! Compares the position of two atoms.
        template< class TTYPE > bool operator < ( const AtomBase<TTYPE> &_atom ) const;
      private:
        //! Serializes an atom.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
          { _ar & pos; _ar & type; }
#       ifdef LADA_WITH_LNS
          //! To load and save to xml-like input.
          template<class T_ARCHIVE> bool lns_access(T_ARCHIVE &_ar, const unsigned int _version);
#       endif
    };

#   ifdef LADA_WITH_LNS
      //! To load and save to xml-like input.
      template<class T_TYPE> template<class T_ARCHIVE>
        bool AtomBase<T_TYPE> lns_access(T_ARCHIVE &_ar, const unsigned int _version) 
        {
          namespace lns = LaDa :: load_n_save;
          return _ar & 
                 (
                   lns::section("Atom")
                     << lns::option( "pos", lns::tag=lns::required, lns::action=pos,
                                     lns::help="Cartesian position in Anstrom." )
                     << lns::option( "type", lns::tag=lns::required, lns::action=type,
                                     lns::help="Atomic specie, or any string." )
                 );
         }
#   endif


    template<class T_TYPE> template< class TTYPE >
      inline bool AtomBase<T_TYPE> :: operator== (const AtomBase<TTYPE> &_atom) const
      {
        return     math::eq( pos[0], _atom.pos[0] ) 
               and math::eq( pos[1], _atom.pos[1] ) 
               and math::eq( pos[2], _atom.pos[2] ); 
      }

    template<class T_TYPE> 
      inline bool AtomBase<T_TYPE> :: equal (const AtomBase<t_Type> &_atom) const
        {
          return     math::eq(_atom.pos[0], pos[0] )
                 and math::eq(_atom.pos[1], pos[1] )
                 and math::eq(_atom.pos[2], pos[2] )
                 and type == _atom.type;
        };

    //! Prints an atom to a stream.
    template<class T_TYPE>
      std::ostream& operator<<(std::ostream &_stream, AtomBase<T_TYPE> const &_atom)
      {
        return stream << std::fixed << std::setprecision(5);
                      << std::setw(8) << _atom.pos[0] << "  "
                      << std::setw(8) << _atom.pos[1] << "  " 
                      << std::setw(8) << _atom.pos[2] << " -- " 
                      << std::setw(16) << _atom.type; 
      }


    template<class T_TYPE> template< class TTYPE >
      bool AtomBase<T_TYPE> :: operator < ( const AtomBase<TTYPE> &_atom ) const
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

  } // namespace Crystal
} // namespace LaDa
  
#endif
