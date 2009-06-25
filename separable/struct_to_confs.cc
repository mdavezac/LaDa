//
//  Version: $Id: ce.cc 1059 2009-04-12 18:48:23Z Mayeul $
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <struct_to_confs.h>

namespace LaDa
{
  namespace CE
  {
    namespace configurations
    {
      void Bases::Origin :: begin( Bases::t_Structure const& _str )
      {
        structure_ = &_str;
        iterator_ = _str.atoms.begin();
      }
      void Bases::Origin :: end( Bases::t_Structure const& _str )
      {
        structure_ = &_str;
        iterator_ = _str.atoms.end();
      }
      void Bases::Xcoord :: create_neighbor_lists_( Bases::Origin const& _origin )
      {
        //! First finds first and second largest distances.
        Bases::t_Structure::t_Atoms::const_iterator i_atom = _origin.structure_->atoms.begin();
        Bases::t_Structure::t_Atoms::const_iterator i_atom_end = _origin.structure_->atoms.end();
        types::t_real min_dist(-1);
        types::t_real second_min_dist(-1);
        atat::rMatrix3d const inv_cell(!_origin.structure_->cell);
        for(; i_atom != i_atom_end; ++i_atom )
        {
          //! Avoid origin.
          if( i_atom == _origin.iterator_ ) continue;
          atat::rVector3d const frac( inv_cell * ( i_atom->pos - *_origin ) );
          atat::rVector3d const centered
          ( 
            frac(0) - std::floor(frac(0) + 0.00001),
            frac(1) - std::floor(frac(1) + 0.00001),
            frac(2) - std::floor(frac(2) + 0.00001)
          );
          atat::rVector3d const cart( _origin.structure_->cell * centered );
          types::t_real const d( atat::norm2(cart) );
          if( d < min_dist or min_dist < 0 )
          {
            if( min_dist >= 0e0 ) second_min_dist = min_dist;
            min_dist = d;
            continue;
          }
          if( d < second_min_dist  or second_min_dist < 0 ) second_min_dist = d;
        }

        //! Now neighbor lists.
        first_.reset( new t_Neighs );
        second_.reset( new t_Neighs );
        size_t fns(0);
        for(; i_atom != i_atom_end; ++i_atom )
        {
          //! Avoid origin.
          bool const is_origin( i_atom == _origin.iterator_ );

          atat::rVector3d const frac( inv_cell * ( i_atom->pos - *_origin ) );
          atat::rVector3d const centered
          ( 
            frac(0) - std::floor(frac(0) + 0.00001),
            frac(1) - std::floor(frac(1) + 0.00001),
            frac(2) - std::floor(frac(2) + 0.00001)
          );
          atat::rVector3d const cart( _origin.structure_->cell * centered );
          for( int i(-1); i < 2; ++i)
            for( int j(-1); j < 2; ++j)
              for( int k(-1); k < 2; ++k)
              {
                if( is_origin and i == 0 and j == 0 and k == 0 ) continue;
                atat::rVector3d const to_add
                ( 
                  cart + _origin.structure_->cell * atat::rVector3d( i, j, k )
                );
                if( Fuzzy::eq( atat::norm2(to_add), min_dist ) )
                {
                  first_->push_back( to_add );
                  ++fns;
                }
                if( fns > 1 ) continue; // no need for second neighbors.
                if( Fuzzy::eq( atat::norm2(to_add), second_min_dist ) )
                  second_->push_back(to_add);
              }
        }
        //! If more than one first neighbor, no need for second neighbors.
        if( fns > 1 ) boost::shared_ptr<t_Neighs>().swap(second_);
      };

      void Bases::Xcoord :: begin( Origin const& _origin )
      {
        create_neighbor_lists_( _origin );
        iterator_ = first_->begin();
        val_ = *iterator_ / atat::norm(*iterator_);
      }
      void Bases::Xcoord :: end( Xcoord const &_b )
      {
        __DOASSERT( not first_, "End iterator needs a valid begin iterator.\n" );
        *this = _b;
        iterator_ = _b.first_->end();
      }
      
      void Bases::Ycoord :: create_equiv_ys( Bases::Xcoord const &_x )
      {
        boost::shared_ptr< Bases::Xcoord :: t_Neighs >
           neighs( _x.second_ ? _x.second_ : _x.first_ );
        types::t_real maxx(-1);
        Bases::Xcoord :: t_Neighs :: const_iterator i_pos = neighs->begin();
        Bases::Xcoord :: t_Neighs :: const_iterator i_pos_end = neighs->end();
        for(; i_pos != i_pos_end; ++i_pos )
        {
          if( i_pos == _x.iterator_ ) continue;
          types::t_real const d( (*i_pos) * (*_x) );
          if( d > maxx or maxx < 0e0 ) maxx = d;
        }
        equivs_.reset( new std::list<atat::rVector3d> );
        for(i_pos = neighs->begin(); i_pos != i_pos_end; ++i_pos)
          if( Fuzzy::eq( (*i_pos) * (*_x), maxx ) ) equivs_->push_back( *i_pos );
      } 
      
      void Bases::Ycoord :: begin( Bases::Xcoord const &_x )
      {
        create_equiv_ys( _x );
        xval_ = *_x;
        iterator_ = equivs_->begin();
        val_ = (*iterator_) - ( (*iterator_) * xval_ ) * xval_;
        val_ = val_ / atat::norm(val_); 
      }
      void Bases::Ycoord :: end( Ycoord const &_b )
      {
        *this = _b;
        iterator_ = equivs_->end();
      }

      Bases::const_iterator Bases::const_iterator :: begin_( Bases :: t_Structure const& _str )
      {
        oiterator_.begin(_str);
        oiterator_end_.end(_str);
        N_ = _str.atoms.size();
        xiterator_.begin(oiterator_); xiterator_end_.end(xiterator_);
        yiterator_.begin(xiterator_); yiterator_end_.end(yiterator_);

        val_.origin = *oiterator_;
        val_.x = *xiterator_;
        val_.y = *yiterator_;
        val_.z = val_.x ^ val_.z;
        val_.weight = 1e0 / types::t_real( N_*yiterator_.size()*xiterator_.size() );
        return *this;
      }

      Bases::const_iterator Bases::const_iterator :: end_( Bases :: t_Structure const& _str )
      {
        oiterator_.end(_str);
        oiterator_end_ = oiterator_;
        xiterator_ = xiterator_end_;
        yiterator_ = yiterator_end_;
        return *this;
      };

      Bases::const_iterator Bases::begin() const { return const_iterator().begin_(structure_); }
      Bases::const_iterator Bases::end() const { return const_iterator().end_(structure_); } 

    }  // end of configuration namespace.
  } // end of Crystal namespace.
} // namespace LaDa
