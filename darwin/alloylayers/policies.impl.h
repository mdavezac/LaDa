//
//  Version: $Id$
//

#include <algorithm>
#include <sstream>

namespace GA
{
  namespace AlloyLayers
  {
    template< class T_OBJECT >
      void Translate<T_OBJECT> :: translate( const Crystal::Structure& _structure,
                                             t_Object& _object )
      {
        // Clears current object.
        typename t_Object :: t_Container& container = _object.Container();
        container.clear();

        const atat::rVector3d direction( _structure.cell.get_column(0) );
        typedef typename Crystal::Structure::t_Atoms::const_iterator t_iatom;
        typedef typename t_Object :: t_Container :: iterator t_ivar;
        t_iatom i_atom = _structure.atoms.begin();
        t_iatom i_atom_end = _structure.atoms.end();
        types::t_real last_depth( i_atom->pos * direction );
        types::t_unsigned layer_size(0);
        types::t_real c(0);
        for( ; i_atom != i_atom_end; )
        {
          types::t_real new_depth;
          do
          {
            c += i_atom->type > 0e0 ? 1e0: -1e0;
            ++layer_size;
            ++i_atom;
            new_depth = i_atom->pos * direction;
          }
          while( Fuzzy::eq( last_depth, new_depth ) and i_atom != i_atom_end );

          // New layer has been reached
          container.push_back( c / (types::t_real) layer_size );

          // Reinitializing layer discovery.
          layer_size = 0;
          last_depth = new_depth;
          c = 0e0;
        }
        std::cout << _object << "\n";
      }

    template< class T_OBJECT >
      void Translate<T_OBJECT> :: translate( const t_Object& _object,
                                             Crystal::Structure& _structure ) 
      {
        types::t_unsigned (*ptr_to_rng)(const types::t_unsigned& )
            = &eo::random<types::t_unsigned>;
        const atat::rVector3d direction( _structure.cell.get_column(0) );
        typedef typename Crystal::Structure::t_Atoms::iterator t_iatom;
        typedef typename t_Object :: t_Container :: const_iterator t_ivar;
        t_iatom i_atom = _structure.atoms.begin();
        t_iatom i_atom_end = _structure.atoms.end();
        t_iatom i_atom_begin( i_atom );
        t_ivar i_var = _object.Container().begin();
        __DODEBUGCODE( t_ivar i_var_end = _object.Container().end(); )
        types::t_real last_depth( i_atom->pos * direction );
        types::t_unsigned layer_size(0);
        for( ; i_atom != i_atom_end __DODEBUGCODE( and i_var != i_var_end ); )
        {
          types::t_real new_depth;
          do
          {
            ++layer_size;
            ++i_atom;
            new_depth = i_atom->pos * direction;
          }
          while( Fuzzy::eq( last_depth, new_depth ) and i_atom != i_atom_end );

          // New layer has been reached
          // Creates an array with approximate concentration, shuffles, then
          // assigns to atoms.
          const bool is_positive( *i_var > 0 );
          const bool is_negative(not is_positive);
          const types::t_real nb_layers_real
          (
            is_positive ? ( 1e0 - (*i_var) ) * 0.5e0 * types::t_real( layer_size ):
                          ( 1e0 + (*i_var) ) * 0.5e0 * types::t_real( layer_size )
          );
          const size_t nb_layers_size_t( nb_layers_real );
          std::vector<bool> arrangment(layer_size,  is_positive );
          std::fill
          ( 
            arrangment.begin(),
            arrangment.begin() + nb_layers_size_t,
            is_negative
          );
          std::random_shuffle( arrangment.begin(), arrangment.end(), ptr_to_rng );
          std::vector<bool> :: const_iterator i_bool = arrangment.end(); 
          for(; i_atom_begin != i_atom; ++i_atom_begin, ++i_bool )
            i_atom_begin->type = *i_bool ? 1e0: -1e0;

          // Reinitializing layer discovery.
          layer_size = 0;
          last_depth = new_depth;
          i_atom_begin = i_atom;
          ++i_var;
        }
        __ASSERT( i_atom != i_atom_end, "Structure smaller than object.\n" )
        __ASSERT( i_var != i_var_end, "Object smaller than structure.\n" )
      }


    template< class T_OBJECT >
      void Translate<T_OBJECT> :: translate( const std::string& _string, 
                                             t_Object& _object ) 
      {
        typedef typename t_Object :: t_Container :: iterator t_ivar;
        std::istringstream sstr( _string );
        t_ivar i_var = _object.Container().begin();
        __DODEBUGCODE( t_ivar i_var_end = _object.Container().end(); )
        do
        {
          __ASSERT( i_var == i_var_end, "Object smaller than string.\n" )
          sstr >> *i_var;
          ++i_var;
        }
        while( not sstr.eof() );
        __ASSERT( i_var != i_var_end, "String smaller than object.\n" )
      }

    template< class T_OBJECT >
      void Translate<T_OBJECT> :: translate( const t_Object& _object,
                                             std::string &_string ) 
      {
        typedef typename t_Object :: t_Container :: const_iterator t_ivar;
        std::ostringstream sstr( _string );
        t_ivar i_var = _object.Container().begin();
        t_ivar i_var_end = _object.Container().end();
        for(; i_var != i_var_end; ++i_var ) sstr << (*i_var) << " ";
        _string = sstr.str(); 
      }

    template< class T_INDIVIDUAL, class T_TRANSLATE >
      bool initialize( T_INDIVIDUAL &_indiv,
                       Crystal::Structure& _structure,
                       T_TRANSLATE _translate )
                       
      {
        foreach( Crystal::Structure::t_Atom &atom, _structure.atoms )
          atom.type = eo::rng.flip() ? 1e0: -1e0;
        _translate( _structure, _indiv.Object() );
        return true;
      }

  }
}
