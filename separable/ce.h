//
//  Version: $Id: separable.h 626 2008-06-18 22:05:05Z davezac $
//
#ifndef _SEPARABLE_CE_H_
#define _SEPARABLE_CE_H_

#include <lamarck/atom.h>

//! Should come to replace less sensibly named Ising_CE/VA_CE
namespace CE
{
  //! A boolean basis of one true and one false function.
  class BooleanBasis;
  //! \brief A scalar function which takes a boolean arguments.
  class Boolean
  {
    friend class BooleanBasis;
    public:
      //! Type of the argument.
      typedef bool t_Arg;
      //! Type of the return.
      typedef types::t_real t_Return;
      //! Does not have gradient.
      bool const static has_gradient;
     
      //! Constructor 
      Boolean ( bool _which = false ) : has_gradient(false), which( _which ) {}
      //! Destructor
      ~Separable() {}

      //! evaluates the function over a range of positions.
      t_Return operator()( const t_Arg _bool ) const
        { return  _bool == which ? t_Return(0): t_Return(1); }

    protected:
      //! Decides whether this is a true or false function.
      bool which;
  };

  class BooleanBasis : public::array< Boolean, 2 >
  {
    public:
      //! Type of argument.
      typedef t_Basis :: value_type :: t_Arg t_Arg;
      //! Type of return.
      typedef t_Basis :: value_type :: t_Return t_Return;
      //! Does not have gradient.
      bool const static has_gradient;

    public:
      //! Constructor.
      BooleanBasis() : has_gradient( t_Basis :: value_type :: has_gradient )
        { elemns[0].which = true; elemns[1].which = false; }
      //! Destructor.
      ~BooleanBasis() {}
  }

// //! A separable function which takes a boolean arguments.
// class Separable
// {
//   public:
//     //! Type of the argument.
//     typedef bool t_Arg;
//     //! Type of the return.
//     typedef types::t_real t_Return;
//     //! Does not have gradient.
//     bool const static has_gradient;
//    
//     //! Constructor 
//     Separable   ( const t_Atom & _atom ) 
//               : pos( _atom.pos ), type( _atom.type ),
//                 coef( 1.0 ), has_gradient(false) {}
//     //! Destructor
//     ~Separable() {}
//
//     //! evaluates the function over a range of positions.
//     template< class T_ITERATOR >
//     t_Return operator()( T_ITERATOR _first, T_ITERATOR _last );
//
//   protected:
//     //! type of the coefficients.
//     typedef std::pair< t_Return, t_Return > t_Coefficient;
//     //! Type of the map
//     typedef std::vector< t_Coefficient > t_Coefficients;
//     //! A container of coefficients
//     t_Coefficients coefs;
// };
//
// //! A sum of separable CE functions.
// class SumSeparable : public ::Separable::Base< Separable >  {};
// bool Separable::Compare::operator()( const t_Key &_a, const t_Key &_b ) const
// {
//   if ( not fuzzy::eq( a[0], b[0] ) )
//     return fuzzy::le( a[0], b[0] );
//   if ( not fuzzy::eq( a[1], b[1] ) )
//     return fuzzy::le( a[1], b[1] );
//   return fuzzy::le( a[2], b[2] );
// }
//
// void Separable :: add( atat::rVector3d & _pos )
// {
//   // position is already in basis.
//   if( map.find( _pos ) != map.end() ) return;
//   map[ _pos ] = t_Coefficient( 1,1 );
// }

// template< class T_ITERATOR > 
//   Separable::t_Return Separable :: operator() ( T_ITERATOR _first,
//                                                 T_ITERATOR _last ) const
//   {
//     t_Return result( *_first `);
//     t_Coefficients :: const_iterator i_coef = coefs.begin();
//     for(; _first != _last; ++_first, ++i_coef )
//     {
//       __ASSERT( i_coef != coefs.end(), "Inconsistent size.\n" );
//       result *= *_first ? i_coef->first : i_found->second;
//     }
//     return 
//   }
} // end of CE namespace

#endif //  _SEPARABLE_CE_H_
