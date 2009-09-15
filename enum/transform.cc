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
     void Transform :: init(atat::rMatrix3d const &_left, atat::iVector3d const &_smith)
             throw(boost::exception)
     {
       std::cout << "ini: 0\n";
       namespace bt = boost::tuples;
       std::cout << "ini: 1\n";
       size_t const nsites_ = independents_.size();
       std::cout << "ini: 2 " << independents_.size() << "\n";
       size_t const card_(nsites_*_smith(0)*_smith(1)*_smith(2));
       std::cout << "ini: 3 " << card_ << "\n";
       if( card_ < 2 ) return;
       std::cout << "ini: 4 " << nsites_ << "\n";
#      ifdef LADA_DEBUG
         if( independents_.size() != nsites_ ) 
         {
           std::ostringstream sstr;
           sstr << independents_.size() << " != " << nsites_ << " ";
           BOOST_THROW_EXCEPTION( internal() << error_string(sstr.str()));
         }
#      endif
       std::cout << "ini: 5" << _smith << "\n";

       namespace bt = boost::tuples;
       permutations_.clear();
       std::cout << "ini: 6\n";

       // loops over sites.
       t_Independents :: const_iterator i_ind = independents_.begin();
       atat::rMatrix3d const rotation = _left * op * (!_left);
       for(types::t_int d(0); d < types::t_int(nsites_); ++d, ++i_ind)
       {
         std::cout << "ini: 7\n";
         types::t_int permutated_site( d + i_ind->first );
         if( permutated_site < 0 ) permutated_site += types::t_int(nsites_);
         
         atat::rVector3d const t_nd
           = _left * atat::rVector3d(i_ind->second(0), i_ind->second(1), i_ind->second(3));
         atat::rVector3d g;
         std::cout << "ini: 8\n";
         // loops over first _smith coordinate.
         for(size_t i(0); i < _smith(0); ++i)
         {
           std::cout << "ini: 9\n";
           g(0) = i;
           // loops over second _smith coordinate.
           for(size_t j(0); j < _smith(1); ++j)
           {
             std::cout << "ini: 10\n";
             g(1) = j;
             // loops over third _smith coordinate.
             for(size_t k(0); k < _smith(2); ++k)
             {
               std::cout << "ini: 11\n";
               g(2) = k;
               std::cout << "ini: 11.1\n";
               atat::rVector3d const transformed( rotation * g + t_nd );
               std::cout << "ini: 11.2\n";
#              ifdef LADA_DEBUG
               std::cout << "ini: 11.3\n";
                 for( size_t t(0); t < 3; ++t)
                 {
                   std::cout << "ini: 11.4\n";
                   if( not Fuzzy::is_zero(transformed(t) - std::floor(transformed(t)+0.5) - 0.5) ) 
                   {
                     std::cout << "ini: 11.5\n";
                     throw symmetry_not_of_supercell();
                     std::cout << "ini: 11.6\n";
                   }
                 }
#              endif
               std::cout << "ini: 12\n";
               atat::iVector3d translation
               (
                 types::t_int(std::floor(transformed(0)+0.5)) % _smith(0),
                 types::t_int(std::floor(transformed(1)+0.5)) % _smith(1), 
                 types::t_int(std::floor(transformed(2)+0.5)) % _smith(2)  
               );
               if( translation(0) < 0 ) translation(0) += _smith(0);
               if( translation(1) < 0 ) translation(1) += _smith(1);
               if( translation(2) < 0 ) translation(2) += _smith(2);

               std::cout << "ini: 13\n";
               permutations_.push_back( get_index(permutated_site, translation, _smith, card_) );
               std::cout << "ini: 14\n";
             } // over k
               std::cout << "ini: 15\n";
           } // over j
               std::cout << "ini: 16\n";
         } // over i
               std::cout << "ini: 17\n";
       } // over d
       std::cout << "ini: 18\n";
     }

     Transform :: Transform   (Crystal::SymmetryOperator const &_c, Crystal::Lattice const &_lat)
                              throw(boost::exception)
                            : Crystal::SymmetryOperator(_c)
     {
       bool pure_translation_ = (    Fuzzy::neq(op.x[0][0], 1.0) 
                                  or Fuzzy::neq(op.x[0][1], 0.0) 
                                  or Fuzzy::neq(op.x[0][2], 0.0)
                                  or Fuzzy::neq(op.x[1][0], 0.0) 
                                  or Fuzzy::neq(op.x[1][1], 1.0) 
                                  or Fuzzy::neq(op.x[1][2], 0.0)
                                  or Fuzzy::neq(op.x[2][0], 0.0) 
                                  or Fuzzy::neq(op.x[2][1], 0.0) 
                                  or Fuzzy::neq(op.x[2][2], 1.0) );

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
         std::cout << "start\n";
         if(card_ < 2) return _x;
         std::cout << "start 0\n";
#      ifdef LADA_DEBUG
         std::cout << "start 1 " << permutations_.size() << " " << card_ << "\n";
         if(permutations_.size() != card_)
           BOOST_THROW_EXCEPTION( internal() << error_string("permutations_ size is incorrect.") );
         std::cout << "start 2\n";
         if(_flavorbase.size() != card_) BOOST_THROW_EXCEPTION( argument_error());
         std::cout << "start 3\n";
         if(_x >= _flavorbase.back() * _flavorbase[1]) BOOST_THROW_EXCEPTION( integer_too_large() );
         std::cout << "start 4\n";
#      endif

         std::cout << "here\n";
       t_uint result(0);
       FlavorBase::const_reverse_iterator i_flavor = _flavorbase.rbegin();
       std::vector<size_t> :: const_iterator i_perm = permutations_.begin();
       std::vector<size_t> :: const_iterator const i_perm_end = permutations_.end();
         std::cout << "here 2\n";
       for(;i_perm != i_perm_end; ++i_flavor, ++i_perm)
       {
         std::cout << "here 3\n";
         t_uint const flavor( _x / (*i_flavor) );
         _x %= (*i_flavor);

         std::cout << "here 4: " << _flavorbase.size() << " " << *i_perm << "\n";
         result += flavor * _flavorbase[*i_perm];
       } // c
         std::cout << "here 4\n";
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
