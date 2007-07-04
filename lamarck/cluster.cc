#include <string>
#include <iomanip>
#include <algorithm>

#ifndef __PGI
  #include<ext/functional>
  using __gnu_cxx::compose1;
  using __gnu_cxx::compose2;
#else
  #include<functional>
  using std::compose2;
  using std::compose1;
#endif
using  std::bind2nd;

#include "opt/types.h"

#include "atat/findsym.h" 
#include "atat/xtalutil.h" 
#include "cluster.h"


namespace Ising_CE {

  void Cluster :: apply_symmetry(const atat::rMatrix3d &point_op, const atat::rVector3d &trans) 
  {
    if ( vectors.size() < 1 ) 
      return;
    std::vector<atat::rVector3d> :: iterator i_vec = vectors.begin();
    std::vector<atat::rVector3d> :: iterator i_last = vectors.end();
    for(; i_vec != i_last; ++i_vec)
      *i_vec = point_op * (*i_vec) + trans;

    i_vec = vectors.begin();
    atat::rVector3d translate = *i_vec;
    for(; i_vec != i_last; ++i_vec)
      *i_vec -= translate;
  }

  bool Cluster :: equivalent_mod_cell( Cluster &equiv, const atat::rMatrix3d &inv_cell) 
  {
    if ( vectors.size() != equiv.vectors.size() )
      return false;
    if ( vectors.size() == 0 )
      return true;

    std::vector<atat::rVector3d> :: iterator i_vec = vectors.begin();
    std::vector<atat::rVector3d> :: iterator i_vec_last = vectors.end();
    for (; i_vec != i_vec_last; ++i_vec)
    {
      if ( atat::equivalent_mod_cell( *equiv.vectors.begin(), *i_vec, inv_cell) )
      {
        atat::rVector3d shift = (*i_vec) - *(equiv.vectors.begin());
        
        std::vector<atat::rVector3d> :: iterator is_found;
        
        // search for equiv  vector such that|| (*i_vec-shift) - *i_equiv ||  > zero_tolerance
        Real (*ptr_func)(const atat::FixedVector<Real, 3>&) = atat::norm;
        std::vector<atat::rVector3d> :: iterator i2_vec = vectors.begin();
        for(; i2_vec != i_vec_last; ++i2_vec)
        {
          is_found  = std::find_if( equiv.vectors.begin(), equiv.vectors.end(),
                                    compose1( bind2nd( std::less_equal<Real>(), atat::zero_tolerance ), 
                                              compose1( std::ptr_fun(ptr_func),
                                                        bind2nd( std::minus<atat::rVector3d>(),  (*i2_vec-shift)) )
                                            ) ); 
          if ( is_found == equiv.vectors.end() )
            break;
        }
                                      
        if ( i2_vec == i_vec_last ) // if all match
          return true; 
      }
    }
    return false;
  }

  void Cluster :: print_out (  std::ostream &stream)  const
  {
    stream << std::fixed << std::setprecision(2) << std::setw(6);
    stream << " Cluster: " << eci << std::endl;
    
    if (vectors.size() == 0 )
    {
      stream << " No position" << std::endl;
      return;
    }

    std::vector<atat::rVector3d> :: const_iterator i_vec = vectors.begin();
    std::vector<atat::rVector3d> :: const_iterator i_last = vectors.end();
    
    for ( ; i_vec != i_last; ++i_vec)
      stream << " " << ( *i_vec )(0) 
             << " " << ( *i_vec )(1) 
             << " " << ( *i_vec )(2) 
             << std::endl;
  }

  bool Cluster :: Load( const TiXmlElement &_node )
  {
    const TiXmlElement *child;
    types::t_real d; atat::rVector3d vec;

    _node.Attribute("eci", &eci);
    vectors.clear();
    child = _node.FirstChildElement( "spin" );
    for ( ; child; child=child->NextSiblingElement( "spin" ) )
    {
      d = 1.0 ; child->Attribute("x", &d); vec(0) = d;
      d = 1.0 ; child->Attribute("y", &d); vec(1) = d;
      d = 1.0 ; child->Attribute("z", &d); vec(2) = d;
      vectors.push_back(vec);
    }

    return true;
  }

} // namespace Ising_CE



#ifdef _MPI

namespace mpi
{
  template<>
  bool BroadCast :: serialize<Ising_CE::Cluster>( Ising_CE::Cluster &_c )
  {
    if ( not serialize( _c.eci ) ) return false;
    
    types::t_int n = _c.vectors.size();
    if( not serialize( n ) ) return false;
    if ( stage == COPYING_FROM_HERE )
      _c.vectors.resize(n);
    std::vector<atat::rVector3d> :: iterator i_vec = _c.vectors.begin();
    std::vector<atat::rVector3d> :: iterator i_vec_end = _c.vectors.end();
    for(; i_vec != i_vec_end; ++i_vec )
      if ( not serialize( i_vec->x, i_vec->x+3 ) ) return false;

    return true;
  }
}

#endif
