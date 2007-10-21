//
//  Version: $Id$
//
#ifndef _FITNESS_IMPL_H_
#define _FITNESS_IMPL_H_

namespace Fitness
{
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
  inline Base<T_QUANTITYTRAITS, true> :: operator t_Quantity() const
  { 
    if ( not is_valid )
      throw std::runtime_error( " Invalid Fitness !!\n" );
    return quantity;
  }
#ifdef _MPI
  template<class T_QUANTITYTRAITS>
  inline bool Base<T_QUANTITYTRAITS, true> :: broadcast( mpi::BroadCast &_bc )
  {
    return     _bc.serialize( is_valid )
           and t_QuantityTraits::broadcast( quantity, _bc );
  }
#endif

  template<class T_QUANTITYTRAITS>
  inline std::ostream & operator<<( std::ostream &_os, const Base<T_QUANTITYTRAITS, true> &_fit ) 
    {  return _os << (const typename T_QUANTITYTRAITS::t_Quantity&) _fit; }
  template<class T_QUANTITYTRAITS>
  inline std::istream & operator>>( std::istream &_is, Base<T_QUANTITYTRAITS, true> &_fit ) 
    { typename T_QUANTITYTRAITS::t_Quantity d; _is >> d; _fit = d; return _is; } 






  template<class T_QUANTITYTRAITS>
  bool Base<T_QUANTITYTRAITS, false> :: operator<( const Base & _f) const
  {
    if( quantity.size() != _f.quantity.size () )
      throw std::runtime_error("quantities of different size in Multi-Objective Fitness!?\n");
    typename t_Quantity :: const_iterator i_scalar1 = quantity.begin();
    typename t_Quantity :: const_iterator i_scalar1_end = quantity.end();
    typename t_Quantity :: const_iterator i_scalar2 = _f.quantity.begin();
    for(; i_scalar1 != i_scalar1_end; ++i_scalar1, ++i_scalar2 )
      if (  t_ScalarQuantityTraits :: less( *i_scalar2, *i_scalar1 ) ) return false;
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
      if (  t_ScalarQuantityTraits :: greater( *i_scalar2, *i_scalar1 ) ) return false;
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
      if (  t_ScalarQuantityTraits :: equal( *i_scalar2, *i_scalar1 ) ) return false;
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
  inline Base<T_QUANTITYTRAITS, false> :: operator const t_Quantity& () const
  { 
    if ( not is_valid )
      throw std::runtime_error( " Invalid Fitness !!\n" );
    return quantity;
  }
#ifdef _MPI
  template<class T_QUANTITYTRAITS>
  inline bool Base<T_QUANTITYTRAITS, false> :: broadcast( mpi::BroadCast &_bc )
  {
    return     _bc.serialize( is_valid )
           and t_QuantityTraits::broadcast( quantity, _bc )
           and Base<T_QUANTITYTRAITS, false>::t_ScalarFitness::broadcast( _bc );
  }
#endif

  template<class T_QUANTITYTRAITS>
  inline std::ostream & operator<<( std::ostream &_os, const Base<T_QUANTITYTRAITS, false> &_fit ) 
  {
    typename T_QUANTITYTRAITS::t_Quantity :: const_iterator i_q = _fit.quantity.begin();
    typename T_QUANTITYTRAITS::t_Quantity :: const_iterator i_q_end = _fit.quantity.end();
    _os << _fit.quantity.size() << " ";
    for(; i_q != i_q_end; ++i_q )
      _os << *i_q << " ";
    return _os << ( const typename Base<T_QUANTITYTRAITS,false>::t_ScalarFitness& ) _fit;
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
    return _is >> ( typename Base<T_QUANTITYTRAITS, false>::t_ScalarFitness& ) _fit;
  }

}
#endif
