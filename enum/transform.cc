//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <crystal/compare_sites.h>
#include <crystal/smith.h>

#include "exceptions.h"
#include "transform.h"


namespace LaDa
{

  namespace enumeration
  {
     void Transform :: init( Crystal::t_SmithTransform const &_transform ) throw(boost::exception)
     {
       namespace bt = boost::tuples;
       atat::rMatrix3d const &left = bt::get<0>(_transform);
       atat::iVector3d const &smith = bt::get<1>(_transform);
       size_t const nsites_ = independents_.size();
       size_t const card_(nsites_*smith(0)*smith(1)*smith(2));

       namespace bt = boost::tuples;
       permutations_.clear();

       // loops over sites.
       t_Independents :: const_iterator i_ind = independents_.begin();
       atat::rMatrix3d const rotation = left * op * (!left);
       for(types::t_int d(0); d < types::t_int(nsites_); ++d, ++i_ind)
       {
         types::t_int permutated_site( d + i_ind->first );
         if( permutated_site < 0 ) permutated_site += types::t_int(nsites_);
         
         atat::rVector3d const t_nd
           = left * atat::rVector3d(i_ind->second(0), i_ind->second(1), i_ind->second(3));
         atat::rVector3d g;
         // loops over first smith coordinate.
         for(size_t i(0); i < smith(0); ++i)
         {
           g(0) = i;
           // loops over second smith coordinate.
           for(size_t j(0); j < smith(0); ++j)
           {
             g(1) = j;
             // loops over third smith coordinate.
             for(size_t k(0); k < smith(0); ++k)
             {
               g(2) = k;
               atat::rVector3d const transformed( rotation * g + t_nd );
#              ifdef LADA_DEBUG
                 for( size_t t(0); t < 3; ++t)
                   if( not Fuzzy::is_zero(transformed(t) - std::floor(transformed(t)+0.5) - 0.5) ) 
                     BOOST_THROW_EXCEPTION( symmetry_not_of_supercell() );
#              endif
               atat::iVector3d translation
               (
                 types::t_int(std::floor(transformed(0)+0.5)) % smith(0),
                 types::t_int(std::floor(transformed(1)+0.5)) % smith(1), 
                 types::t_int(std::floor(transformed(2)+0.5)) % smith(2)  
               );
               if( translation(0) < 0 ) translation(0) += smith(0);
               if( translation(1) < 0 ) translation(1) += smith(1);
               if( translation(2) < 0 ) translation(2) += smith(2);

               permutations_.push_back( get_index(permutated_site, translation, smith, card_) );
             } // over k
           } // over j
         } // over i
       } // over d
     }

     Transform :: Transform   (Crystal::SymmetryOperator const &_c, Crystal::Lattice const &_lat)
                              throw(boost::exception)
                            : Crystal::SymmetryOperator(_c)
     {
       bool pure_translation_ = 
         (    Fuzzy::neq(op.x[0][0], 1.0) or Fuzzy::neq(op.x[0][1], 0.0) or Fuzzy::neq(op.x[0][2], 0.0)
           or Fuzzy::neq(op.x[1][0], 0.0) or Fuzzy::neq(op.x[1][1], 1.0) or Fuzzy::neq(op.x[1][2], 0.0)
           or Fuzzy::neq(op.x[2][0], 0.0) or Fuzzy::neq(op.x[2][1], 0.0) or Fuzzy::neq(op.x[2][2], 1.0) );

       if( pure_translation_ )
       {
         independents_.resize( _lat.sites.size(), t_Independent(0, _c.trans) );
         return;
       }
       atat::rMatrix3d const inv_cell(!_lat.cell);
       std::vector<Crystal::Lattice::t_Site> sites; sites.reserve( _lat.sites.size() );
       { // construct list of centered sites.
         Crystal::Lattice::t_Sites::const_iterator i_site = _lat.sites.begin();
         Crystal::Lattice::t_Sites::const_iterator const i_site_end = _lat.sites.end();
         for(; i_site != i_site_end; ++i_site)
         {
           atat::rVector3d const frac( inv_cell * i_site->pos );
           atat::rVector3d const centered
           (
             frac(0) - std::floor( frac(0) + 0.5 ),
             frac(1) - std::floor( frac(1) + 0.5 ),
             frac(2) - std::floor( frac(2) + 0.5 ) 
           );
           sites.push_back(*i_site); sites.back().pos = centered;
         }
       } 

       { // finds d_{N,d} and t_{N,d}
         std::vector<Crystal::Lattice::t_Site> :: const_iterator i_site_begin = sites.begin();
         std::vector<Crystal::Lattice::t_Site> :: const_iterator i_site = i_site_begin;
         std::vector<Crystal::Lattice::t_Site> :: const_iterator const i_site_end = sites.end();
         atat::rMatrix3d const rotation(inv_cell * op * _lat.cell);
         for(size_t i(0); i_site != i_site_end; ++i_site, ++i)
         {
           atat::rVector3d const frac( rotation * i_site->pos );
           atat::rVector3d const centered
           (
             frac(0) - std::floor( frac(0) + 0.5 ),
             frac(1) - std::floor( frac(1) + 0.5 ),
             frac(2) - std::floor( frac(2) + 0.5 ) 
           );
           Crystal::CompareSites compsites( *i_site );
           compsites.pos = centered;
           std::vector<Crystal::Lattice::t_Site> :: const_iterator i_found 
             = std::find_if( i_site_begin, i_site_end, compsites );
           if( i_found == i_site_end ) 
             BOOST_THROW_EXCEPTION
             (
                symmetry_not_of_lattice()
                  << error_string("Could not find equivalent position.") 
             );
           if( not compsites(i_found->type) )
             BOOST_THROW_EXCEPTION
             ( 
                symmetry_not_of_lattice()
                  << error_string("Equivalent positions do not have equivalent occupations.") 
             );

           // found equivalent sites at this point. 
           // computes d_{N,d} as an integer which sends a d integer values to its transform.
           types::t_int const d_nd = (i_found - i_site_begin) - types::t_int(i);
           // computes translation vector t_{N,d} (in the original, non-centered lattice).
           atat::rVector3d const t_nd( _lat.cell * frac - (_lat.sites[d_nd+i].pos-_lat.sites[i].pos) );
           // pushes into 
           independents_.push_back( t_Independent(d_nd, t_nd + SymmetryOperator::trans) );
         }
       }
     }

     t_uint Transform::operator()(t_uint _x, FlavorBase const &_flavorbase) const
     {
#      ifdef LADA_DEBUG
         if( card_ < 2 ) return _x;
         if( permutations_.size() != card_ )
           BOOST_THROW_EXCEPTION( internal() << error_string("permutations_ size is incorrect.") );
         if( _flavorbase.size() != card_ )
           BOOST_THROW_EXCEPTION( internal() << error_string("_flavorbase size is incorrect.") );
         if( _x >= _flavorbase.back() * _flavorbase[1] )
           BOOST_THROW_EXCEPTION( internal() << error_string("Argument _x is out of range.") );
#      endif

       t_uint result(0);
       FlavorBase::const_reverse_iterator i_flavor = _flavorbase.rbegin();
       std::vector<size_t> :: const_iterator i_perm = permutations_.begin();
       std::vector<size_t> :: const_iterator const i_perm_end = permutations_.begin();
       for(;i_perm != i_perm_end; ++i_flavor)
       {
         t_uint const flavor( _x / (*i_flavor) );
         _x %= (*i_flavor);

         result += flavor * _flavorbase[*i_perm];
       } // c
       return result;
     }

    boost::shared_ptr< std::vector<Transform> > create_transforms( Crystal::Lattice const &_lat )
      {
        using namespace Crystal;
        boost::shared_ptr< std::vector<SymmetryOperator> > symops( get_symmetries(_lat) );
        boost::shared_ptr< std::vector<Transform> > result( new std::vector<Transform> );
        result->reserve( symops->size() );
        foreach(SymmetryOperator const &symop, *symops)
          result->push_back(Transform(symop, _lat));
        return result;
      }
  } // namespace Crystal
} // namespace LaDa
