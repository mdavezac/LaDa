//
//  Version: $Id$
//
#ifndef _BITSTRING_H_
#define _BITSTRING_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<string>
#include<sstream>
#include <functional>

#include<tinyxml/tinyxml.h>

#include<opt/types.h>
#include <mpi/mpi_object.h>


/** \ingroup Genetic
 * @{ */
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
      ~Object() {};
    
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
      const_iterator begin() const
        { return bitstring.begin();  }
      //! Return a constant iterator to the last character + 1 in the string
      const_iterator end() const
        { return bitstring.end();  }
      //! Return an iterator to the first character in the string
      iterator begin() 
        { return bitstring.begin();  }
      //! Return an iterator to the last character + 1 in the string
      iterator end()
        { return bitstring.end();  }
      
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

#ifdef _MPI
      /** \ingroup MPI
       *  \brief Serializes the bitstring of BitString::Object. */
       bool broadcast ( mpi::BroadCast &_bc )
       {
         return _bc.serialize_container( bitstring );
       }
#endif
  };

  //! \brief Dumps a BitString::Object to a stream.
  //! \details This routine is generally used to do print-out to such outlets as
  //           Print::xmg or Print::out. It is \b not used for XML
  //           input/output. Hence it is not used for saving or restarting. You
  //           can overide this routine with little danger.
  template<class T_CONT>
    std::ostream& operator<<(std::ostream &_stream, const Object<T_CONT> &_o);
  //! \brief Dumps a BitString::Object to an std::string instance.
  //! \details This routine is used for XML input/output. Override at your own risk. 
  //! \note The implementation is a complete repeat of BitString::operator<<(
  //        std::ostream&, const BitString::Object<T_CONT>&). It does not use
  //        BitString::operator<<( std::ostream&, const
  //        BitString::Object<T_CONT>&) either directly or inderectly. There
  //        can be no colateral damage to overriding BitString::operator<<(
  //        std::ostream&, const BitString::Object<T_CONT>& ).
  template<class T_CONT>
    void operator<<(std::string &_str, const Object<T_CONT> &_o);
  //! \brief reads a bitstring obect from an std::string.
  //! \details This routine is a mirror of BitString::operator<<(std::string &,
  //           const Object<T_CONT> &). Override at your own risk.
  template<class T_CONT>
    void operator<<(Object<T_CONT> &_o, const std::string &_str );

  //! Returns true if the spin is up, eg \a _t \> 0
  template< class T_TYPE > inline bool spin_up( T_TYPE _t ) { return _t > T_TYPE(0); }
  //! Returns true if spin is up, \a _t is true
  template<> inline bool spin_up<bool>(bool _t ) { return _t; }
  //! Returns true if the spin is down, eg !(\a _t \> 0)
  template< class T_TYPE > inline bool spin_down( T_TYPE _t ) 
    { return  not spin_up<T_TYPE>(_t); }
  //! Returns true if the spin is down, eg \a _t is false
  template<> inline bool spin_down<bool>( bool _t ) { return not _t; }

  //! Flips the spin from up to down, or down to up
  template< class T_TYPE >
    void inline flip( T_TYPE &_t ) { _t = -_t; }
  //! Flips the spin from up to down, or down to up. bool flavor.
  template<>
    void inline flip< bool >( bool &_t ) { _t = not _t; }

  template< class T_CONTAINER >
    inline void Object<T_CONTAINER> :: mask( types::t_unsigned _start,
                                             types::t_unsigned _end)
    {
      if ( _end > bitstring.size() ) _end = bitstring.size();
      std::for_each( bitstring.begin() + _start, bitstring.begin() + _end,
                     std::ptr_fun( &flip< typename T_CONTAINER::value_type > ) );  
    }
    
  //! \brief Implements a bitstring crossover operator.
  //! \details Crossover creates a chilf bitstring from two parents bitstrings.
  //!          For each character in the bitstring, a (possibly biased) coin is
  //!          flipped.  If the coin returns true, then the child inherit that
  //!          character by from one parent. If it returns false, the child
  //!          inherit that character from the other parent. The coin toss is
  //!          repeated independantly for each character.
  //! \xmlinput 
  //!   \code
  //        <Crossover rate="0.5"/>
  //!   \endcode
  //!   The attribute \a rate gives the bias of the coin toss. \a rate="0.5"
  //!   implies no bias. If the attribute does not exist, then rate is set to
  //!   0.5 by default.
  template< class T_OBJECT >
    class Crossover 
    {
      public:
        //! Type of the object. Presumably BitString::Object derived...
        typedef T_OBJECT t_Object;
      public:
        //! \brief The bias of the coin toss
        //! \details Crossover::rate = 0.5 implies that there is no bias. The child will
        //!          inherit from each parent with similar probability. For
        //!          Crossover::rate \< 0.5, the child is more likely to inherit from
        //!          parent \a _parent in
        //!          Crossover::operator()(t_Object &_off, const t_Object& _parent).
        //!          If Crossover::rate \> 0.5, the child is more likely to inherit from
        //!          parent \a _offspring in
        //!          Crossover::operator()( t_Object &_off, t_Object &_parent). 
        //! \pre Crossover::rate must be in the intervall (0, 1). Endpoints will
        //!      result in cloning one or the other parent.
        types::t_real rate;

      public:
        //! Constructor
        Crossover() : rate (0.5) {}
        //! Constructor, initializes BitString::rate from \a _rate.
        Crossover( types::t_real &_rate ) : rate ( _rate ) {}
        //! Constructor and Initializer. Loads Crossover::rate from XML.
        Crossover( const TiXmlElement &_node ) : rate ( 0.5 ) { Load( _node ); }
        //! Copy Constructor
        Crossover( const Crossover<t_Object> &_k ) : rate (_k.rate) {}
        //! Destructor
        ~Crossover() {}

        //! Loads the value of Crossover::rate from XML, as an attribute of \a _node.
        bool Load( const TiXmlElement &_node )
        {
          _node.Attribute("rate", &rate);
          if ( rate <= types::tolerance ) rate = 0.5;
          if ( 1.0 - rate <= types::tolerance ) rate = 1.0;
          return true;
        }

        //! \brief Binary crossover operator.
        //! \details Creates one single child through a standard bitstring
        //!          crossover. See the description of Crossover and
        //!         Crossover::rate for details.
        //! \param[in,out] _off First parent on input, child on output.
        //! \param[in] _parent Second parent on input.
        bool operator()( t_Object &_off, const t_Object &_parent );
        //! \brief Quadratic crossover operator.
        //! \details Creates two children through a standard bitstring
        //!          crossover of two parents. The first child is the one
        //!          described in BitString::Crossover and
        //!          BitString::Crossover::rate. The second child is its
        //!          "complement". In other words, where the first child
        //!          inherits a character from one parent, the second child
        //!          inherits that same character from the other parent. 
        //! \param[in,out] _off1 First parent on input, first child on output.
        //! \param[in,out] _off2 Second parent on input, second child on output.
        bool operator()( t_Object &_off1, t_Object &_off2 );
        //! \brief Returns "Crossover rate = ? " where ? is replaced by the
        //!        value of Crossover::rate.
        std::string print_out() const
        {
          std::ostringstream sstr;
          sstr << "Crossover rate = " << rate;
          return sstr.str();
        }
        //! \brief Returns "BitString::Crossover"
        std::string className() const { return "BitString::Crossover"; }
    };

  template< class T_OBJECT >
    inline bool Crossover<T_OBJECT> :: operator()( t_Object &_off, const t_Object &_parent )
    {
      if( rate <= types::tolerance  ) return false;
      if( 1.0 - rate <= types::tolerance )
      {
        _off = _parent;
        return true;
      }

      typename t_Object :: iterator  i_var_off = _off.bitstring.begin();
      typename t_Object :: iterator  i_var_off_end = _off.bitstring.end();
      typename t_Object :: const_iterator  i_var_par = _parent.bitstring.begin();
      typename t_Object :: const_iterator  i_var_par_end = _parent.bitstring.end();
      bool has_changed = false;
      for(; i_var_off != i_var_off_end and i_var_par != i_var_par_end; ++i_var_off, ++i_var_par )
      {
        if( not rng.flip( rate ) ) continue;
        *i_var_off = *i_var_par;
        has_changed = true;
      }
      return has_changed;
    }
  template< class T_OBJECT >
    inline bool Crossover<T_OBJECT> :: operator()( t_Object &_off1, t_Object &_off2 )
    {
      if( rate <= types::tolerance )
      {
        _off2 = _off1;
        return true;
      }
      if( 1.0 - rate <= types::tolerance )
      {
        _off1 = _off2;
        return true;
      }

      typename t_Object :: iterator  i_var_1 = _off1.bitstring.begin();
      typename t_Object :: iterator  i_var_1_end = _off1.bitstring.end();
      typename t_Object :: iterator  i_var_2 = _off2.bitstring.begin();
      typename t_Object :: iterator  i_var_2_end = _off2.bitstring.end();
      for(; i_var_1 != i_var_1_end and i_var_2 != i_var_2_end; ++i_var_1, ++i_var_2 )
        if( rng.flip( rate ) ) *i_var_1 = *i_var_2;
        else                   *i_var_2 = *i_var_1;
      return true;
    }
   
  //! \brief Implements a bitstring mutation operator.
  //! \details A coin, possibly (preferably!) biased, is tossed for each
  //           character of string. If the toss returns true, the character is
  //           flipped ( see Bitstring::flip() ). If the toss returns false,
  //           the character stays the same.
  //! \xmlinput 
  //!   \code
  //        <Mutation rate="0.1"/>
  //!   \endcode
  //!   The attribute \a rate gives the bias of the coin toss. \a rate="0.1"
  //!   implies no bias. If the attribute does not exist, then rate is set to
  //!   0.1 by default.
  template< class T_OBJECT >
    class Mutation 
    {
      public:
        //! Type of the object. Presumably BitString::Object derived...
        typedef T_OBJECT t_Object; 
      public:
        //! \brief The bias of the coin toss
        //! \details Mutation::rate = 0.5 implies that there is no bias. Each character
        //           has 50\% chance of being flipped. For Mutation::rate \<
        //           0.5, each character is less likely to be flipped. For
        //           Mutation::rate > 0.5, each character is more likely to be
        //          flipped.
        //! \pre Mutation::rate must be in the intervall (0, 1). Endpoints will
        //!      result in no change and complete inversion respectively.
        types::t_real rate;

      public:
        //! Constructor
        Mutation() : rate (0.1) {}
        //! Constructor, initializes BitString::rate from \a _rate.
        Mutation( types::t_real &_rate ) : rate ( _rate ) {}
        //! Constructor and Initializer. Loads Crossover::rate from XML.
        Mutation( const TiXmlElement &_node ) : rate ( 0.1 ) { Load(_node); }
        //! Copy Constructor
        Mutation( const Mutation<t_Object> &_k ) : rate (_k.rate) {}
        //! Destructor
        ~Mutation() {}
          
        //! Loads the value of Crossover::rate from XML, as an attribute of \a _node.
        bool Load(const TiXmlElement &_node )
        {
          _node.Attribute("rate", &rate);
          if ( rate <= types::tolerance ) rate = 0.1;
          if ( 1.0 - rate <= types::tolerance ) rate = 1.0;
          return true;
        }

        //! \brief Mutation operator. 
        bool operator()( t_Object &_o);
        //! \brief Returns "Mutation rate = ? " where ? is replaced by the
        //!        value of Mutation::rate.
        std::string print_out() const
        {
          std::ostringstream sstr;
          sstr << "Mutation rate = " << rate;
          return sstr.str();
        }
        //! \brief Returns "BitString::Mutation"
        std::string className() const { return "BitString::Mutation"; }
    };


  template< class T_OBJECT >
    inline bool Mutation<T_OBJECT> :: operator()( t_Object &_o )
    {
      if( rate <= types::tolerance  ) return false;
      if( 1.0 - rate <= types::tolerance ) rate = 1.0;

      typename t_Object :: iterator  i_var = _o.bitstring.begin();
      typename t_Object :: iterator  i_var_end = _o.bitstring.end();
      bool has_changed = false;
      for(; i_var != i_var_end; ++i_var )
      {
        if( not rng.flip( rate ) ) continue;
        flip< typename t_Object :: t_Container :: value_type >( *i_var );
        has_changed = true;
      }
      return has_changed;
    }

  //! Compares two bitsrings lexicographically
  template<class T_INDIVIDUAL>
  class Bitorder
  {
    public:
      //! functor
      bool operator()(const T_INDIVIDUAL &_i1, const T_INDIVIDUAL &_i2)
      {
        typedef typename T_INDIVIDUAL :: t_IndivTraits :: t_Object :: t_Container 
                                      :: const_iterator t_iterator;
        t_iterator i_var1 = _i1.Object().begin();
        t_iterator i_var1_end = _i1.Object().end();
        t_iterator i_var2 = _i2.Object().begin();
        for(; i_var1 != i_var1_end; ++i_var1, ++i_var2 )
          if( *i_var1 < *i_var2 ) return true;
          else if( *i_var1 > *i_var2 ) return false;
        return false;
      }
  };
}
/* @} */

#include "bitstring.impl.h"

#endif
