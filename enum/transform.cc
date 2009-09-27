//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <crystal/compare_sites.h>
#include <crystal/smith.h>
#include <atat/is_int.h>

#include "exceptions.h"
#include "transform.h"


namespace LaDa
{

  namespace enumeration
  {
     void Transform :: init(atat::rMatrix3d const &_left, atat::iVector3d const &_smith)
     {
       namespace bt = boost::tuples;
       nsites_ = independents_.size();
       card_ = nsites_*_smith(0)*_smith(1)*_smith(2);
       if( card_ < 2 ) return;
#      ifdef LADA_DEBUG
         if( independents_.size() != nsites_ ) 
         {
           std::ostringstream sstr;
           sstr << independents_.size() << " != " << nsites_ << " ";
           BOOST_THROW_EXCEPTION( internal() << error_string(sstr.str()));
         }
#      endif

       namespace bt = boost::tuples;
       permutations_.clear();

       // loops over sites.
       t_Independents :: const_iterator i_ind = independents_.begin();
       atat::rMatrix3d const rotation = _left * op * (!_left);
       bool non_trivial = false;
       for(types::t_int d(0), u(card_-1); d < types::t_int(nsites_); ++d, ++i_ind)
       {
         types::t_int const permutated_site( i_ind->first );
         
         atat::rVector3d const t_nd = _left * i_ind->second;
         atat::rVector3d g;
         // loops over first _smith coordinate.
         for(size_t i(0); i < _smith(0); ++i)
         {
           g(0) = i;
           // loops over second _smith coordinate.
           for(size_t j(0); j < _smith(1); ++j)
           {
             g(1) = j;
             // loops over third _smith coordinate.
             for(size_t k(0); k < _smith(2); ++k, --u)
             {
               g(2) = k;
               atat::rVector3d const transformed( rotation * g + t_nd );
#              ifdef LADA_DEBUG
                 if( not atat::is_integer(transformed) )
                 {
                   throw symmetry_not_of_supercell();
                   BOOST_THROW_EXCEPTION( symmetry_not_of_supercell() );
                 }
#              endif
               atat::iVector3d translation
               (
                 types::t_int(std::floor(transformed(0)+0.001)) % _smith(0),
                 types::t_int(std::floor(transformed(1)+0.001)) % _smith(1), 
                 types::t_int(std::floor(transformed(2)+0.001)) % _smith(2)  
               );
               if( translation(0) < 0 ) translation(0) += _smith(0);
               if( translation(1) < 0 ) translation(1) += _smith(1);
               if( translation(2) < 0 ) translation(2) += _smith(2);

               size_t const index(get_index(permutated_site, translation, _smith, card_));
               permutations_.push_back(index);
               non_trivial |= (u!=index);
             } // over k
           } // over j
         } // over i
       } // over d
       is_trivial_ = not non_trivial;
     }

     Transform :: Transform   (Crystal::SymmetryOperator const &_c, Crystal::Lattice const &_lat)
                            : Crystal::SymmetryOperator(_c)
     {
       bool pure_translation_ = (     Fuzzy::eq(op.x[0][0], 1.0) 
                                  and Fuzzy::eq(op.x[0][1], 0.0) 
                                  and Fuzzy::eq(op.x[0][2], 0.0)
                                  and Fuzzy::eq(op.x[1][0], 0.0) 
                                  and Fuzzy::eq(op.x[1][1], 1.0) 
                                  and Fuzzy::eq(op.x[1][2], 0.0)
                                  and Fuzzy::eq(op.x[2][0], 0.0) 
                                  and Fuzzy::eq(op.x[2][1], 0.0) 
                                  and Fuzzy::eq(op.x[2][2], 1.0) );

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
         std::vector<Crystal::Lattice::t_Site> :: const_iterator const i_site_begin = sites.begin();
         std::vector<Crystal::Lattice::t_Site> :: const_iterator i_site = i_site_begin;
         std::vector<Crystal::Lattice::t_Site> :: const_iterator const i_site_end = sites.end();
         types::t_int const nsites( sites.size() );
         atat::rMatrix3d const rotation(inv_cell * op * _lat.cell);
         for(size_t i(0); i_site != i_site_end; ++i_site, ++i)
         {
           atat::rVector3d const frac( rotation * i_site->pos + inv_cell * trans );
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
           types::t_int const d_nd = ( i_found - i_site_begin );
           // computes translation vector t_{N,d} (in the centered lattice).
           atat::rVector3d const
             t_nd( SymmetryOperator::operator()(_lat.cell*i_site->pos) - _lat.cell * centered );

           // pushes into 
           independents_.push_back( t_Independent(d_nd, t_nd) );
         }
       }
     }

     t_uint Transform::operator()(t_uint _x, FlavorBase const &_flavorbase) const
     {
         if(card_ < 2) return _x;
#      ifdef LADA_DEBUG
         if(permutations_.size() != card_)
           BOOST_THROW_EXCEPTION( internal() << error_string("permutations_ size is incorrect.") );
         if(_flavorbase.size() != card_) BOOST_THROW_EXCEPTION( argument_error());
         if(_x >= _flavorbase.back() * _flavorbase[1]) BOOST_THROW_EXCEPTION( integer_too_large() );
#      endif

       t_uint result(0);
       FlavorBase::const_reverse_iterator i_flavor = _flavorbase.rbegin();
       std::vector<size_t> :: const_iterator i_perm = permutations_.begin();
       std::vector<size_t> :: const_iterator const i_perm_end = permutations_.end();
       for(;i_perm != i_perm_end; ++i_flavor, ++i_perm)
       {
         t_uint const flavor( _x / (*i_flavor) );
         _x %= (*i_flavor);

         if(flavor) result += flavor * _flavorbase[*i_perm];
       } // c
       return result;
     }

    boost::shared_ptr< std::vector<Transform> > create_transforms( Crystal::Lattice const &_lat )
      {
        using namespace Crystal;
        boost::shared_ptr< std::vector<SymmetryOperator> >
          symops( get_space_group_symmetries(_lat) );
        boost::shared_ptr< std::vector<Transform> > result( new std::vector<Transform> );
        result->reserve( symops->size() );
        foreach(SymmetryOperator const &symop, *symops)
          result->push_back(Transform(symop, _lat));
        return result;
      }
  } // namespace Crystal
} // namespace LaDa
