//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<list>
#include<boost/lambda/lambda.hpp>

#include <crystal/lattice.h>
#include <crystal/smith.h>
#include <crystal/symmetry_operator.h>

#include "find_all_cells.h"


namespace LaDa
{

  namespace enumeration
  {
    typedef std::vector<atat::rMatrix3d> t_Container;

    void unique_add( t_Container::value_type const &_supercell, 
                     t_Container::value_type const &_cell, 
                     boost::shared_ptr< std::vector<Crystal::SymmetryOperator> > const &_symops,
                     t_Container &_out );

    boost::shared_ptr<t_Container> find_all_cells(Crystal::Lattice const &_lattice, size_t _nmax)
    {
      boost::shared_ptr< t_Container > result(new t_Container);
      t_Container::value_type cell; cell.zero();

      boost::shared_ptr< std::vector<Crystal::SymmetryOperator> > symops
        = Crystal::get_space_group_symmetries(_lattice);

      // Iterates over values of a such that a * b * c == _nmax
      for(size_t a(1); a <= _nmax; ++a)
      {
        if(_nmax % a != 0) continue;
        size_t const Ndiv_a = _nmax/a;
        cell(0,0) = a;
        // Iterates over values of a such that a * b * c == _nmax
        for(size_t b(1); b <= Ndiv_a; ++b) 
        {
          cell(1,1) = b;
          if(Ndiv_a % b != 0) continue;
          // Iterates over values of a such that a * b * c == _nmax
          size_t const c( Ndiv_a/b);
          cell(2,2) = c;
          if( a * b *c != _nmax ) std::cout << a << " " << b << " " << c << "\n";
          for(size_t d(0); d < b; ++d) 
          {
            cell(1,0) = d;
            for(size_t e(0); e < c; ++e) 
            {
              cell(2,0) = e;
              for(size_t f(0); f < c; ++f) 
              {
                cell(2,1) = f;
                unique_add(cell, _lattice.cell, symops, *result);
              } // f
            } // e
          } // d
        } // b
      } // a

      return result;
    }

    struct IsInt
    {
      IsInt(t_Container::value_type const &_op): op(_op) {}
      IsInt(IsInt const &_c) : op(_c.op) {}
      bool operator()(t_Container::value_type const& _mat)
      {
        for(size_t i(0); i < 3; ++i)
          for(size_t j(0); j < 3; ++j)
          {
            types::t_real a(0);
            for(size_t k(0); k < 3; ++k)
              a += op(i,k) * _mat(k,j);
            if( not Fuzzy::is_zero(std::floor(a+0.1) - a) ) return false;
          }
        return true;
      }
      t_Container::value_type const op;
    };

    void unique_add( t_Container::value_type const &_supercell, 
                     t_Container::value_type const &_cell, 
                     boost::shared_ptr< std::vector<Crystal::SymmetryOperator> > const &_symops,
                     t_Container &_out )
    {
      std::vector<Crystal::SymmetryOperator>::const_iterator i_sym = _symops->begin();
      std::vector<Crystal::SymmetryOperator>::const_iterator const i_sym_end = _symops->end();
      t_Container::value_type const inv(!(_cell*_supercell));
      for(; i_sym != i_sym_end; ++i_sym)
        if( _out.end() != std::find_if(_out.begin(), _out.end(), IsInt(inv*i_sym->op*_cell)) ) 
          break;
      if( i_sym == i_sym_end ) _out.push_back(_supercell);
    }

    boost::shared_ptr< std::vector<SmithGroup> >
      create_smith_groups( Crystal::Lattice const &_lattice,
                           boost::shared_ptr<t_Container> const & _s )
      {
        namespace bt = boost::tuples;
        namespace bl = boost::lambda;
        typedef std::vector<SmithGroup> t_Return;
        boost::shared_ptr<t_Return> result(new t_Return);
        SmithGroup :: Supercell supercell;

        t_Container::const_iterator i_first = _s->begin();
        t_Container::const_iterator const i_end = _s->end();
        for(; i_first != i_end; ++i_first)
        {
          Crystal::t_SmithTransform const
            smith_transform(Crystal::get_smith_transform(_lattice.cell, _lattice.cell*(*i_first)));
          supercell.transform = bt::get<0>(smith_transform);
          supercell.hermite = *i_first;
          t_Return :: iterator i_found
            = std::find_if(result->begin(), result->end(),
                           bl::_1 == bl::constant(bt::get<1>(smith_transform)));
          if( i_found == result->end() )
          {
            result->push_back( SmithGroup(bt::get<1>(smith_transform)) );
            result->back().supercells.push_back(supercell);
          }
          else i_found->supercells.push_back(supercell);
        }
        return result;
      }
    
  } // namespace enumeration

} // namespace LaDa
