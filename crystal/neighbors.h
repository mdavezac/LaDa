//
//  Version: $Id: confsplit.h 844 2008-11-08 01:22:54Z davezac $
//
#ifndef _LADA_CRYSTAL_NEIGHBORS_H_
#define _LADA_CRYSTAL_NEIGHBORS_H_

#include "LaDaConfig.h"

#include <vector>
#include <utility>
#include <algorithm>

#include <boost/lambda/lambda.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <math/fuzzy.h>

#include "structure.h"

namespace LaDa
{

# ifdef LADA_DO_PYTHON
  //! \cond
  namespace Python
  {
    struct Neighbors;
  }
  //! \endcond
# endif

  namespace Crystal
  {
     //! Holds neigbor data.
     struct Neighbor
     {
       Neighbor() {}
       Neighbor( const Neighbor &_c ) : index(_c.index), pos(_c.pos), distance(_c.distance) {}
       //! Index to atom in structure.
       size_t index;
       //! Position with respect to atom.
       math::rVector3d pos;
       //! Distance from atom.
       types::t_real distance;
     };

     inline bool operator<( Neighbor const& _a, Neighbor const& _b )
       { return _a.distance < _b.distance; }

     class Neighbors
     {
#        ifdef LADA_DO_PYTHON
         friend class Python::Neighbors;
#        endif
         //! Type of the list of neighbours.
         typedef std::list<Neighbor> t_Neighbors;
       public:
         //! Iterator over neighbors.
         typedef t_Neighbors::const_iterator const_iterator;
         //! The number of first neighbors to compute.
         size_t nmax;
         //! The origin from which to compute neighbors.
         math::rVector3d origin;

         //! Constructor.
         Neighbors   (size_t _nmax = 0, math::rVector3d const &_vec = math::rVector3d(0,0,0) )
                   : nmax(_nmax), origin(_vec) {};
         //! returns iterator to first neighbor list.
         const_iterator begin() const { return neighbors_.begin(); }
         //! constructs first neighbor list and returns first iterator.
         template<class T_TYPE> 
           const_iterator begin(Crystal::TStructure<T_TYPE> const& _str) 
             { create_neighbors_list_(_str); return begin(); }
         const_iterator begin(Lattice const &_lat)
         {
           Crystal::Structure structure;
           structure.cell = _lat.cell;
           structure.scale = _lat.scale;
           size_t i(0);
           foreach( Lattice::t_Site const &site, _lat.sites )
           {
             Crystal::Structure::t_Atom atom(site.pos, 0);
             atom.site = i;
             structure.atoms.push_back(atom);
             ++i;
           }
           return begin(structure);
         }
         //! returns end of neighbors list.
         const_iterator end() const { return neighbors_.end(); }
         //! Returns size of neighbors list.
         size_t size() const { return neighbors_.size(); }


       private:
         //! Private copy Constructor.
         Neighbors(Neighbors const &_c) {}

         //! Creates list of atoms.
         template<class T_TYPE>
           void create_neighbors_list_(Crystal::TStructure<T_TYPE> const& _str);
         //! List of neighbors.
         t_Neighbors neighbors_;
     };

     template<class T_TYPE>
       void Neighbors :: create_neighbors_list_(Crystal::TStructure<T_TYPE> const& _structure)
       {
         const types::t_int N( _structure.atoms.size() );
         neighbors_.clear();
         
         math::rMatrix3d const inv_cell( !_structure.cell );
         typedef Crystal::TStructure<T_TYPE> t_Structure;
         Neighbor neighbor;
         types::t_real const volume( std::abs(_structure.cell.determinant()) );
         size_t list_max_size(nmax+2);
retry: 
         size_t size(0);
         neighbor.index = 0;
         // Finds out how far to look.
         math::rVector3d const a0( _structure.cell.col(0) );
         math::rVector3d const a1( _structure.cell.col(1) );
         math::rVector3d const a2( _structure.cell.col(2) );
         types::t_real const max_norm
           = std::max( a0.norm(), std::max(a1.norm(), a2.norm()) );
         types::t_real const r
         ( 
           std::pow
           (
             std::max(1e0, types::t_real(list_max_size) / types::t_real(N)),
             0.3333333333333
           )
         );
         types::t_int n0( std::max(1.0, std::ceil(r*max_norm*a1.cross(a2).norm()/volume)) );
         types::t_int n1( std::max(1.0, std::ceil(r*max_norm*a2.cross(a0).norm()/volume)) );
         types::t_int n2( std::max(1.0, std::ceil(r*max_norm*a0.cross(a1).norm()/volume)) );
         while( n0 * n1 * n2 * 8 * N < list_max_size ) { ++n0; ++n1; ++n2; }


         typename t_Structure::t_Atoms::const_iterator i_atom = _structure.atoms.begin();
         typename t_Structure::t_Atoms::const_iterator i_atom_end = _structure.atoms.end();
         for(; i_atom != i_atom_end; ++i_atom, ++neighbor.index ) 
         {
           math::rVector3d const frac( inv_cell * (i_atom->pos - origin) );
           math::rVector3d const centered
           ( 
             frac(0) - std::floor( frac(0) + 0.500000001e0 ),
             frac(1) - std::floor( frac(1) + 0.500000001e0 ),
             frac(2) - std::floor( frac(2) + 0.500000001e0 ) 
           );
           for( types::t_int x(-n0); x <= n0; ++x )
             for( types::t_int y(-n1); y <= n1; ++y )
               for( types::t_int z(-n2); z <= n2; ++z )
               {
                  neighbor.pos = _structure.cell * ( centered + math::rVector3d(x,y,z) );
                  neighbor.distance = neighbor.pos.norm();
                  if( math::is_zero( neighbor.distance ) ) continue;
       
                  t_Neighbors :: iterator i_found 
                  (
                    std::find_if
                    (
                      neighbors_.begin(), neighbors_.end(),
                      boost::lambda::constant(neighbor) < boost::lambda::_1  
                    ) 
                  );
       
                  
                  if( i_found != neighbors_.end() ) 
                  {
                    neighbors_.insert(i_found, neighbor);
       
                    if( size < list_max_size ) ++size;
                    else              neighbors_.pop_back();
                    continue;
                  }
                  
                  if( size == list_max_size ) continue;
       
                  ++size;
                  neighbors_.push_back( neighbor );
               }
         } // loop over atoms.
         
         // Removes atoms beyon nth position which are not at same distance as nth position.
         t_Neighbors::iterator i_last = neighbors_.begin();
         t_Neighbors::iterator const i_end = neighbors_.end();
         LADA_ASSERT( neighbors_.size() > nmax, "Supercell too small.\n");
         size_t i(1);
         for(; i < nmax; ++i, ++i_last );
         LADA_DOASSERT( i_last != i_end, "Supercell too small.\n");
         types::t_real const dist(i_last->distance);
         for(++i_last; i_last != i_end; ++i_last, ++i) 
           if( math::gt(i_last->distance, dist) ) break;
         if( i_last == i_end )
         {
           neighbors_.clear();
           list_max_size += 20;
           goto retry;
         }
         neighbors_.resize(i);
       };

  } // end of Crystal namespace.
} // namespace LaDa

#endif
