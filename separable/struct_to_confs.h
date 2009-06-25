//
//  Version: $Id: confsplit.h 844 2008-11-08 01:22:54Z davezac $
//
#ifndef _SEPARABLES_STRUCT_TO_CONFS_H_
#define _SEPARABLES_STRUCT_TO_CONFS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <utility>

#include <opt/types.h>
#include <opt/debug.h>
#include <crystal/structure.h>

namespace LaDa
{
  namespace CE
  {
    namespace configurations
    {
      struct Basis 
      {
        atat::rVector3d origin;
        atat::rVector3d x;
        atat::rVector3d y;
        atat::rVector3d z;
        types::t_real weight;
      };

      class Bases 
      {
        public:
          //! Iterator over origins.
          class Origin;
          //! Iterator over x coordinates.
          class Xcoord;
          //! Iterator over y coordinates.
          class Ycoord;

          //! Type of structures used for splitting confs.
          typedef Crystal::Structure t_Structure;
          //! Iterator over basis.
          struct const_iterator;

          //! Constructor.
          Bases( t_Structure const &_str) : structure_(_str) {}
          //! Copy Constructor.
          Bases( Bases const &_c) : structure_(_c.structure_) {} 

          //! First basis.
          const_iterator begin() const;
          //! End of basis.
          const_iterator end() const;

        protected:
          t_Structure const &structure_;
      };

      //! Comparison of origin iterators.
      bool operator==( Bases::Origin const&, Bases::Origin const& );

      class Bases::Origin
      {
        friend class Xcoord;
        friend class Ycoord;
        friend bool operator==( Origin const&, Origin const& );
        public:
          //! Return on deref.
          typedef atat::rVector3d const& value_type;
          //! Constructor.
          Origin() {}
          //! Copy Constructor.
          Origin( Origin const& _c ) : iterator_(_c.iterator_) {};

          //! Deref.
          value_type operator*() const { return iterator_->pos; }
          //! Pre-increment operator.
          Origin& operator++() { ++iterator_; return *this; }
          //! Post-increment operator.
          Origin  operator++(int) { ++iterator_; return Origin(*this); }

          //! Points to first origin.
          void begin( t_Structure const& _structure );
          //! Points to end origin.
          void end( t_Structure const& _structure );

        protected:
          //! Pointer to structure. Ugly.
          Bases::t_Structure const * structure_;
          //! Iterator itself.
          Bases::t_Structure::t_Atoms::const_iterator iterator_;
      };
      
      inline bool operator==( Bases::Origin const& _a, Bases::Origin const& _b )
        { return _a.iterator_ != _b.iterator_; }
      inline bool operator!=( Bases::Origin const& _a, Bases::Origin const& _b )
        { return not(_a == _b); }

      //! Comparison of x iterators.
      bool operator==( Bases::Xcoord const&, Bases::Xcoord const& );
      
      
      //! Type of the iterator over x axis.
      class Bases :: Xcoord 
      {
        friend class Ycoord;
        friend bool operator==( Xcoord const&, Xcoord const& );
        //! Type of the vector of neighbors.
        typedef std::vector< atat::rVector3d > t_Neighs;
        public:
          typedef atat::rVector3d const& value_type;
      
          //! Constructor.
          Xcoord() {}
          //! Constructor.
          Xcoord( Origin const& _origin ) { begin(_origin); }
          //! Copy Constructor.
          Xcoord   ( Xcoord const &_c )
                 : first_(_c.first_), second_(_c.second_),
                   iterator_(_c.iterator_), val_(_c.val_) {}
      
          //! Deref.
          value_type operator*() const { return val_; }
          //! Pre-increment operator.
          Xcoord &operator++()
            { ++iterator_; val_ = *iterator_/atat::norm(*iterator_); return *this; }
          //! Post-increment operator.
          Xcoord operator++(int) { return Xcoord( ++(*this) ); }

          //! Sets to begin.
          void begin( Origin const& _origin );
          //! Sets to end.
          void end( Xcoord const &);
          //! Number of equivalent configurations.
          size_t size() const { return first_ ? first_->size(): 0; }

        protected:
          //! Finds first and second neihbors.
          void create_neighbor_lists_( Origin const& _origin );
          //! first neighbors.
          boost::shared_ptr< t_Neighs > first_;
          //! second neighbors, if needed.
          boost::shared_ptr< t_Neighs > second_;
          //! The current iterator.
          t_Neighs::const_iterator iterator_;
          //! Normalized x vector.
          atat::rVector3d val_;
      };

      inline bool operator==( Bases::Xcoord const& _a, Bases::Xcoord const& _b )
        { return _a.iterator_ == _b.iterator_; }
      inline bool operator!=( Bases::Xcoord const& _a, Bases::Xcoord const& _b )
        { return not(_a == _b); }

      //! Comparison of y iterators.
      bool operator==( Bases::Ycoord const& _a, Bases::Ycoord const& _b );

      class Bases :: Ycoord
      {
        friend bool operator==( Ycoord const&, Ycoord const& );
        public:
          typedef atat::rVector3d const& value_type;

          //! Constructor.
          Ycoord() {}
          //! Constructor.
          Ycoord( Xcoord const& _x ) { begin(_x); }
          //! Copy constructors
          Ycoord( Ycoord const& _c ) : equivs_(_c.equivs_), val_(_c.val_), xval_(_c.xval_) {}

          //! Deref.
          value_type operator*() const { return val_; }
          //! Pre-increment operator.
          Ycoord &operator++()
          {
            ++iterator_;
            val_ = (*iterator_) - ( (*iterator_) * xval_ ) * xval_;
            val_ = val_ / atat::norm(val_); 
            return *this; 
          }
          //! Post-increment operator.
          Ycoord operator++(int) { return Ycoord( ++(*this) ); }

          //! Sets to begin.
          void begin( Xcoord const& _x );
          //! Sets to end.
          void end(Ycoord const&);
          //! Number of equivalent configurations.
          size_t size() const { return equivs_ ? first_->equivs_(): 0; }

        protected:
          //! Creates list of equivalent y-positions.
          void create_equiv_ys( Xcoord const &_x );
          //! list of equivalent ys.
          boost::shared_ptr< std::list<atat::rVector3d> > equivs_;
          //! Current iterator.
          std::list<atat::rVector3d> :: const_iterator iterator_;
          //! Holds current dereference value.
          atat::rVector3d val_;
          //! Holds current x coordinate.
          atat::rVector3d xval_;
      };

      inline bool operator==( Bases::Ycoord const& _a, Bases::Ycoord const& _b )
        { return _a.iterator_ == _b.iterator_; }
      inline bool operator!=( Bases::Ycoord const& _a, Bases::Ycoord const& _b )
        { return not(_a == _b); }

      //! Comparison of basis iterators.
      bool operator==( Bases::const_iterator const& _a, Bases::const_iterator const& _b );

      class Bases::const_iterator
      {
        friend class Bases;
        friend bool operator==( const_iterator const&, const_iterator const& );
        public:
          typedef Basis const& value_type;

          //! Constructor
          const_iterator() {};
          //! Copy Constructor.
          const_iterator   ( const_iterator const& _c )
                         : val_(_c.val_),
                           oiterator_(_c.oiterator_), oiterator_end_(_c.oiterator_end_),
                           xiterator_(_c.xiterator_), xiterator_end_(_c.xiterator_end_),
                           yiterator_(_c.yiterator_), yiterator_end_(_c.yiterator_end_) {}

          //! Dereference.
          value_type operator*() const { return val_; }
          //! Pre-increment operator.
          const_iterator &operator++()
          {
            if( oiterator_ == oiterator_end_ ) return *this;

            ++yiterator_;
            if( yiterator_ != yiterator_end_ ) 
            {
              val_.y = *yiterator_;
              val_.z = val_.x ^ val_.z;
              val_.weight = 1e0 / types::t_real( N_*yiterator_->size()*xiterator_->size() );
              return *this;
            }
              
            ++xiterator_;
            if( xiterator_ != xiterator_end_ )
            {
              yiterator_ = Bases::Ycoord( xiterator_ );
              yiterator_.begin(xiterator_); yiterator_end_.end(yiterator_);
              val_.x = *xiterator_;
              val_.y = *yiterator_;
              val_.z = val_.x ^ val_.z;
              val_.weight = 1e0 / types::t_real( N_*yiterator_->size()*xiterator_->size() );
              return *this;
            }
           
            ++oiterator_;
            if( oiterator_ == oiterator_end_ ) return *this;

            xiterator_.begin(oiterator_); xiterator_end_.end(xiterator_);
            yiterator_.begin(xiterator_); yiterator_end_.end(yiterator_);
            val_.origin = *oiterator_;
            val_.x = *xiterator_;
            val_.y = *yiterator_;
            val_.z = val_.x ^ val_.z;
            val_.weight = 1e0 / types::t_real( N_*yiterator_->size()*xiterator_->size() );
            return *this; 
          }
          //! Post-increment operator.
          const_iterator operator++(int) { return const_iterator( ++(*this) ); }

        private:
          //! initializes to first.
          const_iterator begin_( Bases :: t_Structure const& _str );
          //! initializes to first.
          const_iterator end_( Bases :: t_Structure const& _str );
          //! Value to return when dereferenced.
          Basis val_;
          //! origin iterator.
          Bases::Origin oiterator_;
          //! origin iterator end.
          Bases::Origin oiterator_end_;
          //! x iterator.
          Bases::Xcoord xiterator_;
          //! x iterator end. 
          Bases::Xcoord xiterator_end_;
          //! y iterator.
          Bases::Ycoord yiterator_;
          //! y iterator.
          Bases::Ycoord yiterator_end_;
          //! Size of the structure.
          size_t N_;
      };

      inline bool operator==( Bases::const_iterator const& _a, Bases::const_iterator const& _b )
        { return _a.oiterator_ == _b.oiterator_; }
      inline bool operator!=( Bases::const_iterator const& _a, Bases::const_iterator const& _b )
        { return not(_a == _b); }
    }  // end of configuration namespace.
  } // end of Crystal namespace.
} // namespace LaDa

#endif
