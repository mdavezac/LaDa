//
//  Version: $Id$
//
#ifndef _FITNESS_H_
#define _FITNESS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<vector>
#include<math.h>
#include<stdexcept>

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#ifdef _MPI
#include <mpi/mpi_object.h>
#endif
#include <print/xmg.h>
#include <print/manip.h>

#include<opt/traits.h>

namespace Fitness
{
  template<class T_QUANTITYTRAITS,
           bool IS_SCALAR = T_QUANTITYTRAITS::is_scalar >
  class Base {};

  template<class T_QUANTITYTRAITS >
  class Base<T_QUANTITYTRAITS, true>
  {
    typedef Base<T_QUANTITYTRAITS, true> t_This;
    template<class TQUANTITYTRAITS>
      friend inline std::istream & operator>>( std::ostream &_is,
                                               const Base<TQUANTITYTRAITS, true> &_fit );
    template<class TQUANTITYTRAITS>
      friend inline std::ostream & operator<<( std::istream &_is,
                                               Base<TQUANTITYTRAITS, true> &_fit );
    public:
      typedef T_QUANTITYTRAITS t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_Quantity t_Quantity;
      typedef t_This t_ScalarFitness;

    protected:
      t_Quantity quantity;
      bool is_valid;

    public:
      Base() : is_valid( false )  {}
      Base( const Base & _c ) : quantity( _c.quantity ), is_valid( _c.is_valid ) {}
      Base( const types::t_real _fit ) : quantity( _fit ), is_valid( true ) {}
      ~Base() {}


      bool operator<(const Base & _f) const
        { return std::abs(quantity - _f.quantity) > types::tolerance and quantity > _f.quantity; }
      bool operator>(const Base & _f) const
        { return std::abs(quantity - _f.quantity) > types::tolerance and quantity < _f.quantity; }
      bool operator==(const Base & _f) const
        { return std::abs(quantity - _f.quantity) < types::tolerance; }

      bool invalid() const { return not is_valid; }
      void invalidate() { is_valid = false; }
      operator t_Quantity() const
      { 
        if ( not is_valid )
          throw std::runtime_error( " Invalid Fitness !!\n" );
        return quantity;
      }
      bool Load( const TiXmlElement & _node );
      bool Save( TiXmlElement & _node ) const;

#ifdef _MPI
      bool broadcast( mpi::BroadCast &_bc )
      {
        return     _bc.serialize( is_valid )
               and t_QuantityTraits::broadcast( quantity, _bc );
      }
#endif 
  };

  template<class T_QUANTITYTRAITS>
  bool Base<T_QUANTITYTRAITS, true> :: Load( const TiXmlElement & _node ) 
  {
    if ( not _node.Attribute( "fitness" ) )
    {
      is_valid = false;
      return false; 
    }
    double d; _node.Attribute( "fitness", &d );
    quantity = (t_Quantity) d;
    is_valid = true;
    return true;
  }
  template<class T_QUANTITYTRAITS>
  bool Base<T_QUANTITYTRAITS, true> :: Save( TiXmlElement & _node ) const
  {
    if ( not is_valid )
    {
      std::cerr << "Trying to Save invalid fitness" << std::endl;
      return false;
    }
    double d = (double) quantity;
    _node.SetDoubleAttribute("fitness", d);
    return true;
  }

  template<class T_QUANTITYTRAITS>
  inline std::ostream & operator<<( std::ostream &_os, const Base<T_QUANTITYTRAITS, true> &_fit ) 
    {  return _os << (typename T_QUANTITYTRAITS::t_Quantity) _fit; }
  template<class T_QUANTITYTRAITS>
  inline std::istream & operator>>( std::istream &_is, Base<T_QUANTITYTRAITS, true> &_fit ) 
    { typename T_QUANTITYTRAITS::t_Quantity d; _is >> d; _fit = d; return _is; } 



  template<class T_QUANTITYTRAITS >
  class Base<T_QUANTITYTRAITS, false> :
        public Base< typename T_QUANTITYTRAITS::t_ScalarQuantityTraits, true >
  {
    typedef Base<T_QUANTITYTRAITS, false> t_This;
    template<class TQUANTITYTRAITS>
      friend std::istream & operator>>( std::istream &_is,
                                        Base<TQUANTITYTRAITS, false> &_fit );
    template<class TQUANTITYTRAITS>
      friend std::ostream & operator<<( std::ostream &_os,
                                        const Base<TQUANTITYTRAITS, false> &_fit );
    public:
      typedef T_QUANTITYTRAITS t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_Quantity t_Quantity;
      typedef typename t_QuantityTraits :: t_ScalarQuantityTraits t_ScalarQuantityTraits;
      typedef Base<t_ScalarQuantityTraits, true> t_ScalarFitness;

    protected:
      typedef Base<t_ScalarQuantityTraits, true> t_Base;


    protected:
      t_Quantity quantity;
      bool is_valid;

    public:
      Base() : is_valid( false )  {}
      Base( const Base & _c ) : t_Base(_c), quantity( _c.quantity ), is_valid( _c.is_valid ) {}
      Base( const t_Quantity &_fit ) : t_Base(), quantity( _fit ), is_valid( true ) {}
      Base( const t_ScalarFitness _fit ) : t_Base( _fit ), is_valid( false ) {}
      ~Base() {}


      bool operator<(const Base & _f) const;
      bool operator>(const Base & _f) const;
      bool operator==(const Base & _f) const;

      bool invalid() const { return not is_valid; }
      void invalidate() { is_valid = false; }
      operator const t_Quantity& () const
      { 
        if ( not is_valid )
          throw std::runtime_error( " Invalid Fitness !!\n" );
        return quantity;
      }
      bool Load( const TiXmlElement & _node );
      bool Save( TiXmlElement & _node ) const;

#ifdef _MPI
      bool broadcast( mpi::BroadCast &_bc )
      {
        return     _bc.serialize( is_valid )
               and t_QuantityTraits::broadcast( quantity, _bc );
      }
#endif 
  };

  template<class T_QUANTITYTRAITS>
  bool Base<T_QUANTITYTRAITS, false> :: operator<( const Base & _f) const
  {
    if( quantity.size() != _f.quantity.size () )
      throw std::runtime_error("quantities of different size in Multi-Objective Fitness!?\n");
    typename t_Quantity :: const_iterator i_scalar1 = quantity.begin();
    typename t_Quantity :: const_iterator i_scalar1_end = quantity.end();
    typename t_Quantity :: const_iterator i_scalar2 = _f.quantity.begin();
    for(; i_scalar1 != i_scalar1_end; ++i_scalar1, ++i_scalar2 )
      if (     *i_scalar2 > *i_scalar1  
           or std::abs( *i_scalar1 - *i_scalar2 ) < types::tolerance ) return false;
    return true;
  }
  template<class T_QUANTITYTRAITS>
  bool Base<T_QUANTITYTRAITS, false> :: operator>( const Base & _f) const
  {
    if( quantity.size() != _f.quantity.size () )
      throw std::runtime_error("quantities of different size in Multi-Objective Fitness!?\n");
    typename t_Quantity :: const_iterator i_scalar1 = quantity.begin();
    typename t_Quantity :: const_iterator i_scalar1_end = quantity.end();
    typename t_Quantity :: const_iterator i_scalar2 = _f.quantity.begin();
    for(; i_scalar1 != i_scalar1_end; ++i_scalar1, ++i_scalar2 )
      if (     *i_scalar2 < *i_scalar1  
           or std::abs( *i_scalar1 - *i_scalar2 ) < types::tolerance ) return false;
    return true;
  }
  template<class T_QUANTITYTRAITS>
  bool Base<T_QUANTITYTRAITS, false> :: operator==( const Base & _f) const
  {
    if( quantity.size() != _f.quantity.size () )
      throw std::runtime_error("quantities of different size in Multi-Objective Fitness!?\n");
    typename t_Quantity :: const_iterator i_scalar1 = quantity.begin();
    typename t_Quantity :: const_iterator i_scalar1_end = quantity.end();
    typename t_Quantity :: const_iterator i_scalar2 = _f.quantity.begin();
    for(; i_scalar1 != i_scalar1_end; ++i_scalar1, ++i_scalar2 )
      if ( std::abs( *i_scalar1 - *i_scalar2 ) > types::tolerance ) return false;
    return true;
  }
  template<class T_QUANTITYTRAITS>
  bool Base<T_QUANTITYTRAITS, false> :: Load( const TiXmlElement & _node )
  {
    if ( not _node.Attribute( "fitness" ) )
    {
      is_valid = false;
      return false; 
    }
    std::istringstream istr; istr.str( _node.Attribute( "fitness" ) );
    istr >> *this;
    is_valid = true;
    return true;
  }
  template<class T_QUANTITYTRAITS>
  bool Base<T_QUANTITYTRAITS, false> :: Save( TiXmlElement & _node ) const
  {
    if ( not is_valid )
    {
      std::cerr << "Trying to Save invalid fitness" << std::endl;
      return false;
    }
    std::ostringstream ostr ; ostr << *this;
    _node.SetAttribute("fitness", ostr.str().c_str() );
    return true;
  }

  template<class T_QUANTITYTRAITS>
  inline std::ostream & operator<<( std::ostream &_os, const Base<T_QUANTITYTRAITS, false> &_fit ) 
  {
    typename T_QUANTITYTRAITS::t_Quantity :: const_iterator i_q = _fit.quantity.begin();
    typename T_QUANTITYTRAITS::t_Quantity :: const_iterator i_q_end = _fit.quantity.end();
    _os << _fit.quantity.size() << " ";
    for(; i_q != i_q_end; ++i_q )
      _os << *i_q << " ";
    return _os;
  }
  template<class T_QUANTITYTRAITS>
  inline std::istream & operator>>( std::istream &_is, Base<T_QUANTITYTRAITS, false> &_fit ) 
  {
    types::t_unsigned n; _is >> n;
    _fit.quantity.clear();
    while( n and _is.good() )
    {
      typename T_QUANTITYTRAITS::t_ScalarQuantity scalar = 0;
      _is >> scalar;
      _fit.quantity.push_back( scalar ); 
      --n;
    }
    return _is;
  }
}


#endif // _FITNESS_H_
