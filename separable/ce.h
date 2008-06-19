#ifndef _SEPARABLE_CE_H_
#define _SEPARABLE_CE_H_

#include <lamarck/atom.h>

//! Should come to replace less sensibly named Ising_CE/VA_CE
namespace CE
{
  //! A CE separable function.
  class Separable
  {
    //! A key comparison functor.
    struct Compare;
    protected:
      //! Type of the keys.
      typedef atat::rVector3d t_Key;

    public:
      //! Type of the argument.
      typedef Ising_CE::Atom t_Arg;
      //! Type of the return.
      typedef types::t_real t_Return;
     
      //! Constructor 
      function   ( const t_Atom & _atom ) 
               : pos( _atom.pos ), type( _atom.type ), coef( 1.0 );
      //! evaluates the function over a range of positions.
      template< class T_ITERATOR >
      t_Return operator()( T_ITERATOR _first, T_ITERATOR _last );
      //! \brief evaluates the function over a single position
      //! \details Returns 0 if _arg is not in basis set.
      t_Return operator()( t_Arg &_arg )

      //! Adds a new basis point
      void add( t_Key &_pos );
      //! Returns true if the position is present in the basis
      bool is_present( t_Key &_pos ) const;
        { return map.find( _arg.pos ) !=  map.end(); }

    protected:
      //! type of the coefficients.
      typedef std::pair< t_Return, t_Return > t_Coefficient;
      //! Type of the map
      typedef std::map< t_Key, t_Coefficient, Compare > t_Map;
      //! A map to the basis values.
      t_Map coefs;
  };

  struct Separable::Compare
  {
    //! The functor itself.
    bool operator()( const t_Key &_a, const t_Key &_b ) const
    {
      if ( not fuzzy::eq( a[0], b[0] ) )
        return fuzzy::le( a[0], b[0] );
      if ( not fuzzy::eq( a[1], b[1] ) )
        return fuzzy::le( a[1], b[1] );
      return fuzzy::le( a[2], b[2] );
    }
  }

  void Separable :: add( atat::rVector3d & _pos )
  {
    // position is already in basis.
    if( map.find( _pos ) != map.end() ) return;
    map[ _pos ] = t_Coefficient( 1,1 );
  }

  template< class T_ITERATOR > 
    Separable::t_Return Separable :: operator() ( T_ITERATOR _first, T_ITERATOR _last ) const
    {
      t_Return result(1);
      t_Map :: const_iterator i_end( map.end() );
      for(; _first != _last; ++_first )
      {
        t_Map :: const_iterator i_found( map.find( _first->pos ) );
        if( i_found ==  i_end ) continue;
        result *= _first->type < t_Return(0) ?
                    i_found->second->first :
                    i_found->second->second;
      }
      return result;
    }
    Separable::t_Return Separable :: operator() ( const t_Arg &_arg ) const
    {
      t_Map :: const_iterator i_found( map.find( _arg.pos ) );
      if( i_found ==  map.end() ) return t_Return(1);
      return _arg.type < t_Return(0) ?
                  i_found->second->first :
                  i_found->second->second;
    }

    Separable::t_Return Separable :: operator() ( const t_Arg &_arg ) const
    {
      t_Map :: const_iterator i_found( map.find( _arg.pos ) );
      if( i_found ==  map.end() ) return t_Return(1);
      return _arg.type < t_Return(0) ?
                  i_found->second->first :
                  i_found->second->second;
    }
} // end of CE namespace

#endif //  _SEPARABLE_CE_H_
