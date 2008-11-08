//
//  Version: $Id$
//
#ifndef _FITNESS_IMPL_H_
#define _FITNESS_IMPL_H_

#include <boost/serialization/serialization.hpp>

#include <opt/debug.h>

namespace LaDa
{
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
      __ASSERT( not is_valid, "Trying to save invalid fitness\n" )
      double d = (double) quantity;
      _node.SetDoubleAttribute("fitness", d);
      return true;
    }
    template<class T_QUANTITYTRAITS>
    inline t_Comparison Scalar<T_QUANTITYTRAITS> :: compare( const t_This &_f ) const
    {
      if( this->operator<(_f) ) return WEAKER;
      if( this->operator==(_f) ) return EQUAL;
      return STRONGER;
    }
    template<class T_QUANTITYTRAITS> template< class ARCHIVE >
    void Scalar<T_QUANTITYTRAITS> :: serialize( ARCHIVE & _ar,
                                                const unsigned int _version)
    {
      _ar & quantity;
      _ar & is_valid;
    }



    template<class T_QUANTITYTRAITS>
    t_Comparison Vectorial<T_QUANTITYTRAITS> :: compare( const t_This & _f) const
    {
      __ASSERT( vec_quantity.size() != _f.vec_quantity.size(),
                   "Comparing quantities of different sizes\n"
                << vec_quantity.size() << " "
                << _f.vec_quantity.size() << "\n" )

      typename t_Quantity :: const_iterator i_scalar1 = vec_quantity.begin();
      typename t_Quantity :: const_iterator i_scalar1_end = vec_quantity.end();
      typename t_Quantity :: const_iterator i_scalar2 = _f.vec_quantity.begin();
      t_Comparison result = EQUAL;
      for(; i_scalar1 != i_scalar1_end; ++i_scalar1, ++i_scalar2 )
      {
        if ( t_ScalarQuantityTraits :: eq( *i_scalar2, *i_scalar1 ) )
         continue;
        else if( t_ScalarQuantityTraits :: geq( *i_scalar2, *i_scalar1 ) )
        {
          if( result == STRONGER ) return INDIFFERENT;
          result = WEAKER;
          continue;
        }
        if( result == WEAKER ) return INDIFFERENT;
        result = STRONGER;
      }

      return result;
    }

    template<class T_QUANTITYTRAITS>
    bool Vectorial<T_QUANTITYTRAITS> :: operator<( const t_This & _f) const
    {
      __ASSERT( vec_quantity.size() != _f.vec_quantity.size(),
                   "Comparing quantities of different sizes\n"
                << vec_quantity.size() << " "
                << _f.vec_quantity.size() << "\n" )

      typename t_Quantity :: const_iterator i_scalar1 = vec_quantity.begin();
      typename t_Quantity :: const_iterator i_scalar1_end = vec_quantity.end();
      typename t_Quantity :: const_iterator i_scalar2 = _f.vec_quantity.begin();
      for(; i_scalar1 != i_scalar1_end; ++i_scalar1, ++i_scalar2 )
        if ( not t_ScalarQuantityTraits :: gt( *i_scalar2, *i_scalar1 ) ) return false;
      return true;
    }
    template<class T_QUANTITYTRAITS>
    bool Vectorial<T_QUANTITYTRAITS> :: operator>( const t_This & _f) const
    {
      __ASSERT( vec_quantity.size() != _f.vec_quantity.size(),
                   "Comparing quantities of different sizes\n"
                << vec_quantity.size() << " "
                << _f.vec_quantity.size() << "\n" )

      typename t_Quantity :: const_iterator i_scalar1 = vec_quantity.begin();
      typename t_Quantity :: const_iterator i_scalar1_end = vec_quantity.end();
      typename t_Quantity :: const_iterator i_scalar2 = _f.vec_quantity.begin();
      for(; i_scalar1 != i_scalar1_end; ++i_scalar1, ++i_scalar2 )
        if ( not t_ScalarQuantityTraits :: le( *i_scalar2, *i_scalar1 ) ) return false;
      return true;
    }
    template<class T_QUANTITYTRAITS>
    bool Vectorial<T_QUANTITYTRAITS> :: operator==( const t_This & _f) const
    {
      __ASSERT( vec_quantity.size() != _f.vec_quantity.size(),
                   "Comparing quantities of different sizes\n"
                << vec_quantity.size() << " "
                << _f.vec_quantity.size() << "\n" )

      typename t_Quantity :: const_iterator i_scalar1 = vec_quantity.begin();
      typename t_Quantity :: const_iterator i_scalar1_end = vec_quantity.end();
      typename t_Quantity :: const_iterator i_scalar2 = _f.vec_quantity.begin();
      for(; i_scalar1 != i_scalar1_end; ++i_scalar1, ++i_scalar2 )
        if ( not t_ScalarQuantityTraits :: eq( *i_scalar2, *i_scalar1 ) ) return false;
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
      __ASSERT( not vec_is_valid, "Trying to save invalid fitness\n" )
      
      std::ostringstream ostr ; ostr << *this;
      _node.SetAttribute("fitness", ostr.str().c_str() );
      return true;
    }
    template<class T_QUANTITYTRAITS> template< class ARCHIVE >
    void Vectorial<T_QUANTITYTRAITS> :: serialize( ARCHIVE & _ar,
                                                   const unsigned int _version)
    {
      _ar & boost::serialization::base_object< t_Base >( *this );
      _ar & vec_quantity;
      _ar & vec_is_valid;
    }


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

  } // namespace Fitness
} // namespace LaDa

#endif
