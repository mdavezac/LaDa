//
//  Version: $Id$
//
#ifndef _FITNESS_IMPL_H_
#define _FITNESS_IMPL_H_

namespace Fitness
{
  template<class T_QUANTITYTRAITS>
  bool Scalar<T_QUANTITYTRAITS> :: Load( const TiXmlElement & _node ) 
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
  bool Scalar<T_QUANTITYTRAITS> :: Save( TiXmlElement & _node ) const
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
  inline Scalar<T_QUANTITYTRAITS> :: operator t_Quantity() const
  { 
    if ( not is_valid )
    {
      std::ostringstream sstr;
      sstr << __LINE__ << ", line: " << __LINE__ << "\n"
           << "Tried to access invalid fitness\n";
      throw std::runtime_error( sstr.str() );
    }
    return quantity;
  }
#ifdef _MPI
  template<class T_QUANTITYTRAITS>
  inline bool Scalar<T_QUANTITYTRAITS> :: broadcast( mpi::BroadCast &_bc )
  {
    return     _bc.serialize( is_valid )
           and t_QuantityTraits::broadcast( quantity, _bc );
  }
#endif






  template<class T_QUANTITYTRAITS>
  bool Vectorial<T_QUANTITYTRAITS> :: operator<( const t_This & _f) const
  {
    if( vec_quantity.size() != _f.vec_quantity.size() )
    {
      std::ostringstream sstr;
      sstr << __LINE__ << ", line: " << __LINE__ << "\n"
           << "Comparing quantities of different sizes\n"
           << vec_quantity.size() << " " << _f.vec_quantity.size() << "\n";
      throw std::runtime_error( sstr.str() );
    }
    typename t_Quantity :: const_iterator i_scalar1 = vec_quantity.begin();
    typename t_Quantity :: const_iterator i_scalar1_end = vec_quantity.end();
    typename t_Quantity :: const_iterator i_scalar2 = _f.vec_quantity.begin();
    for(; i_scalar1 != i_scalar1_end; ++i_scalar1, ++i_scalar2 )
      if ( not t_ScalarQuantityTraits :: greater( *i_scalar2, *i_scalar1 ) ) return false;
    return true;
  }
  template<class T_QUANTITYTRAITS>
  bool Vectorial<T_QUANTITYTRAITS> :: operator>( const t_This & _f) const
  {
    if( vec_quantity.size() != _f.vec_quantity.size () )
    {
      std::ostringstream sstr;
      sstr << __LINE__ << ", line: " << __LINE__ << "\n"
           << "Comparing quantities of different sizes\n"
           << vec_quantity.size() << " " << _f.vec_quantity.size() << "\n";
      throw std::runtime_error( sstr.str() );
    }
    typename t_Quantity :: const_iterator i_scalar1 = vec_quantity.begin();
    typename t_Quantity :: const_iterator i_scalar1_end = vec_quantity.end();
    typename t_Quantity :: const_iterator i_scalar2 = _f.vec_quantity.begin();
    for(; i_scalar1 != i_scalar1_end; ++i_scalar1, ++i_scalar2 )
      if ( not t_ScalarQuantityTraits :: less( *i_scalar2, *i_scalar1 ) ) return false;
    return true;
  }
  template<class T_QUANTITYTRAITS>
  bool Vectorial<T_QUANTITYTRAITS> :: operator==( const t_This & _f) const
  {
    if( vec_quantity.size() != _f.vec_quantity.size () )
    {
      std::ostringstream sstr;
      sstr << __LINE__ << ", line: " << __LINE__ << "\n"
           << "Comparing quantities of different sizes\n"
           << vec_quantity.size() << " " << _f.vec_quantity.size() << "\n";
      throw std::runtime_error( sstr.str() );
    }
    typename t_Quantity :: const_iterator i_scalar1 = vec_quantity.begin();
    typename t_Quantity :: const_iterator i_scalar1_end = vec_quantity.end();
    typename t_Quantity :: const_iterator i_scalar2 = _f.vec_quantity.begin();
    for(; i_scalar1 != i_scalar1_end; ++i_scalar1, ++i_scalar2 )
      if ( not t_ScalarQuantityTraits :: equal( *i_scalar2, *i_scalar1 ) ) return false;
    return true;
  }
  template<class T_QUANTITYTRAITS>
  bool Vectorial<T_QUANTITYTRAITS> :: Load( const TiXmlElement & _node )
  {
    if ( not _node.Attribute( "fitness" ) )
    {
      vec_is_valid = false;
      return false; 
    }
    std::istringstream istr; istr.str( _node.Attribute( "fitness" ) );
    istr >> *this;
    vec_is_valid = true;
    return true;
  }
  template<class T_QUANTITYTRAITS>
  bool Vectorial<T_QUANTITYTRAITS> :: Save( TiXmlElement & _node ) const
  {
    if ( not vec_is_valid )
    {
      std::cerr << "Trying to Save invalid fitness" << std::endl;
      return false;
    }
    std::ostringstream ostr ; ostr << *this;
    _node.SetAttribute("fitness", ostr.str().c_str() );
    return true;
  }
  template<class T_QUANTITYTRAITS>
  inline Vectorial<T_QUANTITYTRAITS> :: operator const t_Quantity& () const
  { 
    if ( not vec_is_valid )
    {
      std::ostringstream sstr;
      sstr << __LINE__ << ", line: " << __LINE__ << "\n"
           << "Tried to access invalid fitness\n";
      throw std::runtime_error( sstr.str() );
    }
    return vec_quantity;
  }
#ifdef _MPI
  template<class T_QUANTITYTRAITS>
  inline bool Vectorial<T_QUANTITYTRAITS> :: broadcast( mpi::BroadCast &_bc )
  {
    return     _bc.serialize( vec_is_valid )
           and t_QuantityTraits::broadcast( vec_quantity, _bc )
           and Vectorial<T_QUANTITYTRAITS>::t_ScalarFitness::broadcast( _bc );
  }
#endif

  template<class T_QUANTITYTRAITS>
  inline std::ostream & operator<<( std::ostream &_os,
                                    const Scalar<T_QUANTITYTRAITS> &_fit ) 
    {  return _os << (const typename T_QUANTITYTRAITS::t_Quantity&) _fit; }
  template<class T_QUANTITYTRAITS>
  inline std::istream & operator>>( std::istream &_is, 
                                    Scalar<T_QUANTITYTRAITS> &_fit ) 
    { typename T_QUANTITYTRAITS::t_Quantity d; _is >> d; _fit = d; return _is; } 

  template<class T_QUANTITYTRAITS>
  inline std::ostream & operator<<( std::ostream &_os,
                                    const Vectorial<T_QUANTITYTRAITS> &_fit ) 
  {
    const typename T_QUANTITYTRAITS::t_Quantity & vec_quantity = _fit;
    typename T_QUANTITYTRAITS::t_Quantity :: const_iterator i_q = vec_quantity.begin();
    typename T_QUANTITYTRAITS::t_Quantity :: const_iterator i_q_end = vec_quantity.end();
    _os << vec_quantity.size() << " ";
    for(; i_q != i_q_end; ++i_q )
      _os << *i_q << " ";
    typedef typename Vectorial<T_QUANTITYTRAITS>::t_ScalarFitness t_Base;
    return _os << (const t_Base&) _fit;
  }
  template<class T_QUANTITYTRAITS>
  inline std::istream & operator>>( std::istream &_is,
                                    Vectorial<T_QUANTITYTRAITS> &_fit ) 
  {
    types::t_unsigned n; _is >> n;
    _fit.clear();
    while( n and _is.good() )
    {
      typename T_QUANTITYTRAITS::t_ScalarQuantity scalar = 0;
      _is >> scalar;
      _fit.push_back( scalar ); 
      --n;
    }
    typedef typename Vectorial<T_QUANTITYTRAITS>::t_ScalarFitness t_Base;
    return _is >> ( t_Base& ) _fit;
  }

}


#endif
