//
//  Version: $Id$
//
#ifndef _LADA_GA_BITSTRING_OBJECT_H_
#define _LADA_GA_BITSTRING_OBJECT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<iostream>
#include<string>
#include<functional>

#include<opt/types.h>

namespace LaDa
{
  namespace GA
  {

    //! \brief Basic components for bitstring optimization.
    //! \details Since lattice-site occupations can easily be mapped onto a
    //!          bitstring, the following classes represent the basis of most, if
    //!          not all, of the decoration searches implemented here. In BitString
    //!          are provided stub classes for objects and bitstring operations.
    //!          Still, it aims to be quite general. As such physics-specific
    //!          components, such as Concentration and Fourier Transform functors
    //!          are left for further derived classes. 
    //
    //!          More specifically, BitString provides an object class with a
    //!          bitstring container, a standard crossover functor, and a standard
    //!          mutation operator.
    //! \see For decoration search derived from this namespace, see GroundState,
    //!          BandGap, and Molecularity
    namespace BitString
    {

      //! \brief Stub class for a bitstring container
      //! \details Object contains a bit-string container and defines a few general
      //!          routines over it. For the sake of generality, this container can
      //!          be anything with random access iterators. It can be truly  a
      //!          collection of binary bits (bool), or it could be a collection of
      //!          types::t_real, as in done most everywhere in LaDa. In theory,
      //!          the "characters" of this bitstring could be a more complex type,
      //!          though what behaviors would be required of such a character, I
      //!          do not know.
      template< class T_CONTAINER = std::vector< types::t_real > >
        struct Object
        {
            friend class boost::serialization::access;
          public:
            //! The container type
            typedef T_CONTAINER t_Container; 
            //! The bit type
            typedef typename t_Container :: value_type t_Type; 
            //! Iterator to the container type
            typedef typename t_Container :: iterator iterator; 
            //! Constant iterator to the container
            typedef typename t_Container :: const_iterator const_iterator;
  
          public:
            //! The bitstring itself
            t_Container bitstring;
  
          public:
            //! \brief Constrcutor
            //! \details  Sizes the bitstring to zero
            Object()  { bitstring.resize(0); }
            //! Copy Constructor
            Object(const Object &_c) : bitstring(_c.bitstring) {};
            //! Constructor. Initializes the bitstring from \a _c
            Object(const t_Container &_c) : bitstring(_c) {};
            //! Destructor
            virtual ~Object() {};
          
            //! Returns true if the bitstrings are equal
            bool operator==( const Object &_c ) const
            {
              return std::equal( bitstring.begin(), bitstring.end(), 
                                  _c.bitstring.begin() ); 
            }
  
            //! Inverts the bitstring in the range [ \a _start, \a _end )
            //! \see BitString::flip()
            void mask( types::t_unsigned _start, types::t_unsigned _end);
            //! Returns the size of the bitstring
            types::t_unsigned size() { return bitstring.size(); }
            //! Returns a reference to the bitstring container
            t_Container& Container() { return bitstring; }
            //! Returns a constant reference to the bitstring container
            const t_Container& Container() const { return bitstring; }
  
            //! Return a constant iterator to the first character in the string
            const_iterator begin() const  { return bitstring.begin();  }
            //! Return a constant iterator to the last character + 1 in the string
            const_iterator end() const  { return bitstring.end();  }
            //! Return an iterator to the first character in the string
            iterator begin()  { return bitstring.begin();  }
            //! Return an iterator to the last character + 1 in the string
            iterator end() { return bitstring.end();  }
            
            //! \brief Returns the "concentration", or average character of the string,
            //! \details With \e N the number of characters in the string and
            //          \f$c_i\f$ the i(th) character of the string, this function
            //          returns \f[ \frac{1}{N} \sum_i c_i \f]. Not that the
            //          characters are cast into types::t_real before summation and
            //          division.
            types::t_real get_concentration() const
            {
              typename t_Container :: const_iterator i_var = bitstring.begin();
              typename t_Container :: const_iterator i_var_end = bitstring.end();
              types::t_real result = 0.0;
              for(; i_var != i_var_end; ++i_var )
                result += *i_var > 0 ? 1.0: -1.0;
              result /= static_cast<types::t_real>(bitstring.size());
              return result;
            }
  
          private:
            //! Serializes an atom.
            template<class ARCHIVE> void serialize( ARCHIVE & _ar,
                                                    const unsigned int _version)
              { _ar & bitstring; }
        };

      //! \brief Dumps a BitString::Object to a stream.
      //! \details This routine is generally used to do print-out to such outlets as
      //!          Print::xmg or Print::out. It is \b not used for XML
      //!          input/output. Hence it is not used for saving or restarting. You
      //!          can overide this routine with little danger.
      template<class T_CONT>
        std::ostream& operator<<(std::ostream &_stream, const Object<T_CONT> &_o)
        {
          typedef typename Object<T_CONT> ::t_Type t_Type;
          typedef typename Object<T_CONT> ::t_Container :: const_iterator const_iterator;
          const_iterator i_var = _o.bitstring.begin();
          const_iterator i_end = _o.bitstring.end();
          for(; i_var != i_end; ++i_var )
              _stream << *i_var << " ";
          return _stream;
        }

      template< class T_CONTAINER >
        inline void Object<T_CONTAINER> :: mask( types::t_unsigned _start,
                                                 types::t_unsigned _end)
        {
          if ( _end > bitstring.size() ) _end = bitstring.size();
          std::for_each( bitstring.begin() + _start, bitstring.begin() + _end,
                         std::negate< typename T_CONTAINER::value_type >() );
        }
    } // namespace BitString.

  } // namespace GA


} // namespace LaDa

#endif
