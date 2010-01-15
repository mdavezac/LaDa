//
//  Version: $Id$
//


namespace LaDa
{
  namespace atomic_potential
  {
    template<class T_STRUCTURE>
      void Bases<T_STRUCTURE>::Origin :: begin( T_STRUCTURE const& _str )
      {
        structure_ = &_str;
        iterator_ = _str.atoms.begin();
        index_ = 0;
      }
    template<class T_STRUCTURE>
      void Bases<T_STRUCTURE>::Origin :: end( T_STRUCTURE const& _str )
      {
        structure_ = &_str;
        iterator_ = _str.atoms.end();
        index_ = _str.atoms.size();
      }
    template<class T_STRUCTURE>
      void Bases<T_STRUCTURE>::Xcoord
              ::create_neighbor_lists_( Bases<T_STRUCTURE>::Origin const& _origin )
      {
        // We only need first neighbors, if there are more than one,
        // and second neighbors if there is only one first neighbor.
        // 13 nearest neighbors should be enough to cover these situations.
        Crystal::Neighbors neighs( 13, *_origin );
        Crystal::Neighbors::const_iterator i_first = neighs.begin( *_origin.structure_ );
        Crystal::Neighbors::const_iterator const i_end = neighs.end();
        types::t_real const min_dist( i_first->distance );
        first_.reset( new t_Neighs );
        for(; i_first != i_end; ++i_first )
        {
          if( math::gt( i_first->distance, min_dist ) ) break;
          
          first_->push_back( i_first->pos );
        }
     
        // if more than one first neighbor, this is all we need.
        __DOASSERT( i_first == i_end, "Could not find second-largest distance.\n" );
        if( first_->size() > 1 ) 
        {
          second_ = first_;
          return;
        }
        // otherwise get second nearest neighbors.
        second_.reset( new t_Neighs );
        types::t_real const second_min_dist( i_first->distance );
        for(; i_first != i_end; ++i_first )
        {
          if( math::gt( i_first->distance, second_min_dist ) ) break;
          
          second_->push_back( i_first->pos );
        }
        __DOASSERT( i_first == i_end, "Could not find third-largest distance.\n" );
      };
     
    template<class T_STRUCTURE>
      void Bases<T_STRUCTURE>::Xcoord :: begin( Origin const& _origin )
      {
        create_neighbor_lists_( _origin );
        iterator_ = first_->begin();
        val_ = iterator_->normalized();
      }
    template<class T_STRUCTURE>
      void Bases<T_STRUCTURE>::Xcoord :: end( Xcoord const &_b )
      {
        __DOASSERT( not _b.first_, "End iterator needs a valid begin iterator.\n" );
        *this = _b;
        iterator_ = _b.first_->end();
      }
      
    template<class T_STRUCTURE>
      void Bases<T_STRUCTURE>::Ycoord :: create_equiv_ys( Bases<T_STRUCTURE>::Xcoord const &_x )
      {
        boost::shared_ptr< typename Bases<T_STRUCTURE>::Xcoord :: t_Neighs >
           neighs( _x.second_ ? _x.second_ : _x.first_ );
        types::t_real maxx(-1);
        typename Bases<T_STRUCTURE>::Xcoord :: t_Neighs :: const_iterator i_pos = neighs->begin();
        typename Bases<T_STRUCTURE>::Xcoord :: t_Neighs :: const_iterator i_pos_end = neighs->end();
        for(; i_pos != i_pos_end; ++i_pos )
        {
          if( i_pos == _x.iterator_ ) continue;
          types::t_real const d( (*i_pos) * (*_x) );
          if( d > maxx or maxx < 0e0 ) maxx = d;
        }
        equivs_.reset( new std::list<Eigen::Vector3d> );
        for(i_pos = neighs->begin(); i_pos != i_pos_end; ++i_pos)
          if( math::eq( (*i_pos) * (*_x), maxx ) ) equivs_->push_back( *i_pos );
      } 
      
    template<class T_STRUCTURE>
      void Bases<T_STRUCTURE>::Ycoord :: begin( Bases<T_STRUCTURE>::Xcoord const &_x )
      {
        create_equiv_ys( _x );
        xval_ = *_x;
        iterator_ = equivs_->begin();
        val_ = ( (*iterator_) - ( (*iterator_) * xval_ ) * xval_ ).normalized();
      }
    template<class T_STRUCTURE>
      void Bases<T_STRUCTURE>::Ycoord :: end( Ycoord const &_b )
      {
        *this = _b;
        iterator_ = equivs_->end();
      }
     
    template<class T_STRUCTURE> typename Bases<T_STRUCTURE>::const_iterator 
      Bases<T_STRUCTURE>::const_iterator::begin_( T_STRUCTURE const& _str )
      {
        oiterator_.begin(_str);
        oiterator_end_.end(_str);
        N_ = _str.atoms.size();
        xiterator_.begin(oiterator_); xiterator_end_.end(xiterator_);
        yiterator_.begin(xiterator_); yiterator_end_.end(yiterator_);
     
        val_.origin = *oiterator_;
        val_.index = oiterator_.index();
        val_.x = *xiterator_;
        val_.y = *yiterator_;
        val_.z = val_.x ^ val_.y;
        val_.weight = 1e0 / types::t_real( N_*yiterator_.size()*xiterator_.size() );
        return *this;
      }
     
    template<class T_STRUCTURE> typename Bases<T_STRUCTURE>::const_iterator 
      Bases<T_STRUCTURE>::const_iterator :: end_( T_STRUCTURE const& _str )
      {
        oiterator_.end(_str);
        oiterator_end_ = oiterator_;
        xiterator_ = xiterator_end_;
        yiterator_ = yiterator_end_;
        return *this;
      };
     
    template<class T_STRUCTURE> 
      typename Bases<T_STRUCTURE>::const_iterator Bases<T_STRUCTURE>::begin() const
        { return const_iterator().begin_(structure_); }
    template<class T_STRUCTURE> 
      typename Bases<T_STRUCTURE>::const_iterator Bases<T_STRUCTURE>::end() const
        { return const_iterator().end_(structure_); } 
     
    template<class T_STRUCTURE> 
      typename Bases<T_STRUCTURE>::const_iterator &Bases<T_STRUCTURE>::const_iterator::operator++()
      {
        if( oiterator_ == oiterator_end_ ) return *this;
     
        ++yiterator_;
        if( yiterator_ != yiterator_end_ ) 
        {
          val_.y = *yiterator_;
          val_.z = val_.x ^ val_.y;
          size_t const u(N_*yiterator_.size()*xiterator_.size());
          val_.weight = 1e0 / types::t_real( u );
          return *this;
        }
          
        ++xiterator_;
        if( xiterator_ != xiterator_end_ )
        {
          yiterator_ = typename Bases<T_STRUCTURE>::Ycoord( xiterator_ );
          yiterator_.begin(xiterator_); yiterator_end_.end(yiterator_);
          val_.x = *xiterator_;
          val_.y = *yiterator_;
          val_.z = val_.x ^ val_.y;
          val_.weight = 1e0 / types::t_real( N_*yiterator_.size()*xiterator_.size() );
          return *this;
        }
       
        ++oiterator_;
        if( oiterator_ == oiterator_end_ ) return *this;
     
        xiterator_.begin(oiterator_); xiterator_end_.end(xiterator_);
        yiterator_.begin(xiterator_); yiterator_end_.end(yiterator_);
        val_.origin = *oiterator_;
        val_.index = oiterator_.index();
        val_.x = *xiterator_;
        val_.y = *yiterator_;
        val_.z = val_.x ^ val_.y;
        val_.weight = 1e0 / types::t_real( N_*yiterator_.size()*xiterator_.size() );
        return *this; 
      }
  } // end of atomic_potential namespace.
} // namespace LaDa
