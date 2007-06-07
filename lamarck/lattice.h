#ifndef _LATTICE_H_
#define _LATTICE_H_


#include <vector>
#include <fstream>
#include <string>
#include <complex>
#include <exception>


#include <tinyxml/tinyxml.h>
#include <opt/types.h>

#include "atat/findsym.h"
#include "atat/vectmac.h"
#include "atat/machdep.h"

#include "atom.h"


namespace Ising_CE {

  struct Lattice
  {
    typedef Atom_Type< std::vector<std::string> > t_Site;
    typedef std::vector<t_Site> t_Sites;
    atat::rMatrix3d cell;
    t_Sites sites;
    atat::SpaceGroup space_group;
    types::t_real scale;

    Lattice() {};
    Lattice( const TiXmlElement &_element )
      { Load( _element ); };
    ~Lattice () {};
    bool Load( const TiXmlElement &_element );
    void find_space_group();
    types::t_unsigned get_nb_sites() const
      { return sites.size(); } 
    types::t_unsigned get_nb_types( types::t_unsigned i ) const
      { return sites[i].type.size(); } 
    types::t_int get_atom_site_index( const atat::rVector3d &_at ) const;
    types::t_int get_atom_site_index( const std::string &_at ) const;
    types::t_int get_atom_type_index( const Ising_CE::Atom &_at ) const;
    types::t_int get_atom_type_index( const std::string &_at ) const;
    std::string get_atom_string( const Ising_CE::Atom &_at ) const
    {
      types::t_int i = get_atom_site_index( _at );
      if ( i == -1 ) return "error";
      if ( get_nb_types(i) == 1 ) return sites[i].type[0];
      return ( std::abs( _at.type - 1.0 ) < atat::zero_tolerance ) ? 
          sites[i].type[0] : sites[i].type[1];  
    }
    const std::string& get_atom_string( const unsigned _s, const unsigned _i ) const
    {
      return sites[_s].type[_i];
    }
    bool convert_StrAtom_to_Atom( const Ising_CE::StrAtom &_in,
                                  Ising_CE::Atom &_out ) const;
    bool convert_Atom_to_StrAtom( const Ising_CE::Atom &_in,
                                  Ising_CE::StrAtom &_out ) const;

    protected:
      types::t_real convert_type_index_to_real( const types::t_unsigned _i ) const
        { return ( _i ) ? 1.0 : -1.0; }
      types::t_unsigned convert_real_to_type_index( const types::t_real _r ) const
        { return ( std::abs(_r + 1.0) < atat::zero_tolerance ) ? 0 : 1; }
    public:
      types::t_real convert_type_index_to_real( const types::t_unsigned _s,
                                                const types::t_unsigned _i ) const
      {
        if ( sites[_s].type.size() == 1 ) 
          return -1;
        return convert_type_index_to_real( _i );
      }
      types::t_unsigned convert_real_to_type_index( const types::t_unsigned _s,
                                                    const types::t_real _r ) const
      {
        if ( sites[_s].type.size() == 1 ) 
          return 0;
        return convert_real_to_type_index( _r );
      }
      void print_out (std::ostream &stream) const;
  };

} // namespace Ising_CE
#endif
