//
//  Version: $Id$
//
#ifndef _LADA_ATOMIC_POTENTIALS_BASES_H_
#define _LADA_ATOMIC_POTENTIALS_BASES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <utility>

#include <opt/types.h>
#include <opt/debug.h>
#include <crystal/structure.h>
#include <crystal/neighbors.h>

namespace LaDa
{
  namespace atomic_potential
  {
    struct Basis 
    {
      //! Index of the origin in structure for which a representation is built.
      size_t index;
      //! Cartesian coordinates of the origin.
      math::rVector3d origin;
      //! Cartesian coordinates of the abscissa.
      math::rVector3d x;
      //! Cartesian coordinates of the ordinate.
      math::rVector3d y;
      //! Cartesian coordinates of the cote.
      math::rVector3d z;
      //! Weight of the basis in the representation.
      types::t_real weight;
    };
    inline std::ostream& operator<<( std::ostream &_stream, Basis const &_basis )
    {
      return _stream << "Origin: " << _basis.origin << " -- "
                     << "weight: " << _basis.weight << "\n"
                     << "x: " << _basis.x << "\n"
                     << "y: " << _basis.y << "\n"
                     << "z: " << _basis.z << "\n";
    }

    template<class T_STRUCTURE>
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
          typedef T_STRUCTURE t_Structure;
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

    template<class T_STRUCTURE>
      class Bases<T_STRUCTURE>::Origin
      {
        friend class Xcoord;
        friend class Ycoord;
        public:
          //! Return on deref.
          typedef math::rVector3d const& value_type;
          //! Return on deref.
          typedef math::rVector3d const* pointer_type;
          //! Constructor.
          Origin() {}
          //! Copy Constructor.
          Origin   ( Origin const& _c )
                 : structure_(_c.structure_), iterator_(_c.iterator_), index_(_c.index_) {};
  
          //! Deref.
          value_type operator*() const { return iterator_->pos; }
          //! Deref.
          pointer_type operator->() const { return &(iterator_->pos); }
          //! Returns the index of the origin.
          size_t index() const { return index_; }
          //! Pre-increment operator.
          Origin& operator++() { ++iterator_; ++index_; return *this; }
          //! Post-increment operator.
          Origin  operator++(int) { return operator++(); }
  
          //! Returns false if origins at same position.
          bool operator!=( Origin const& _b ) const { return iterator_ != _b.iterator_; }
          //! Returns true if origins at same position.
          bool operator==( Origin const& _b ) const { return iterator_ == _b.iterator_; }
          //! Points to first origin.
          void begin( t_Structure const& _structure );
          //! Points to end origin.
          void end( t_Structure const& _structure );
  
        protected:
          //! Pointer to structure. Ugly.
          typename Bases<T_STRUCTURE>::t_Structure const * structure_;
          //! Iterator itself.
          typename Bases<T_STRUCTURE>::t_Structure::t_Atoms::const_iterator iterator_;
          //! index of the origin.
          size_t index_;
      };
    
    //! Comparison of x iterators.
    template<class T_STRUCTURE>
      bool operator==( typename Bases<T_STRUCTURE>::Xcoord const&,
                       typename Bases<T_STRUCTURE>::Xcoord const& );
    
    
    //! Type of the iterator over x axis.
    template<class T_STRUCTURE>
      class Bases<T_STRUCTURE> :: Xcoord 
      {
        friend class Ycoord;
        //! Type of the vector of neighbors.
        typedef std::vector< math::rVector3d > t_Neighs;
        public:
          typedef math::rVector3d const& value_type;
          typedef math::rVector3d const* pointer_type;
      
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
          //! Deref.
          pointer_type operator->() const { return &val_; }
          //! Pre-increment operator.
          Xcoord &operator++()
            { ++iterator_; val_ = iterator_->normalized(); return *this; }
          //! Post-increment operator.
          Xcoord operator++(int) { return Xcoord( ++(*this) ); }
  
          //! Sets to begin.
          void begin( Origin const& _origin );
          //! Sets to end.
          void end( Xcoord const &);
          //! Number of equivalent configurations.
          size_t size() const { return first_ ? first_->size(): 0; }

          //! Returns true if x coordinate at same position.
          bool operator==( Xcoord const& _b ) const { return iterator_ == _b.iterator_; }
          //! Returns false if x coordinate at same position.
          bool operator!=( Xcoord const& _b ) const { return iterator_ != _b.iterator_; }
  
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
          math::rVector3d val_;
      };

    template<class T_STRUCTURE>
      class Bases<T_STRUCTURE> :: Ycoord
      {
        public:
          typedef math::rVector3d const& value_type;
          typedef math::rVector3d const* pointer_type;
      
          //! Constructor.
          Ycoord() {}
          //! Constructor.
          Ycoord( Xcoord const& _x ) { begin(_x); }
          //! Copy constructors
          Ycoord   ( Ycoord const& _c )
                 : equivs_(_c.equivs_), iterator_(_c.iterator_),
                   val_(_c.val_), xval_(_c.xval_) {}
      
          //! Deref.
          value_type operator*() const { return val_; }
          //! Deref.
          pointer_type operator->() const { return &val_; }
          //! Pre-increment operator.
          Ycoord &operator++()
          {
            ++iterator_;
            val_ = ( (*iterator_) - ( (*iterator_) * xval_ ) * xval_ ).normalized();
            return *this; 
          }
          //! Post-increment operator.
          Ycoord operator++(int) { return Ycoord( ++(*this) ); }
      
          //! Sets to begin.
          void begin( Xcoord const& _x );
          //! Sets to end.
          void end(Ycoord const&);
          //! Number of equivalent configurations.
          size_t size() const { return equivs_ ? equivs_->size(): 0; }

          //! Returns true if iterators are at same position.
          bool operator==(  Ycoord const& _b ) const 
            { return iterator_ == _b.iterator_; }
          //! Returns false if iterators are at same position.
          bool operator!=(  Ycoord const& _b ) const
            { return iterator_ != _b.iterator_; }
      
        protected:
          //! Creates list of equivalent y-positions.
          void create_equiv_ys( Xcoord const &_x );
          //! list of equivalent ys.
          boost::shared_ptr< std::list<math::rVector3d> > equivs_;
          //! Current iterator.
          std::list<math::rVector3d> :: const_iterator iterator_;
          //! Holds current dereference value.
          math::rVector3d val_;
          //! Holds current x coordinate.
          math::rVector3d xval_;
      };
      
    template<class T_STRUCTURE>
      class Bases<T_STRUCTURE>::const_iterator
      {
        friend class Bases<T_STRUCTURE>;
        public:
          typedef Basis const value_type;
          typedef value_type& reference;
          typedef value_type* pointer_type;
      
          //! Constructor
          const_iterator() {};
          //! Copy Constructor.
          const_iterator   ( const_iterator const& _c )
                         : val_(_c.val_),
                           oiterator_(_c.oiterator_), oiterator_end_(_c.oiterator_end_),
                           xiterator_(_c.xiterator_), xiterator_end_(_c.xiterator_end_),
                           yiterator_(_c.yiterator_), yiterator_end_(_c.yiterator_end_),
                           N_(_c.N_) {}
      
          //! Dereference.
          reference operator*() const { return val_; }
          //! Deref.
          pointer_type operator->() const { return &val_; }
          //! Pre-increment operator.
          const_iterator &operator++();
          //! Post-increment operator.
          const_iterator operator++(int) { return const_iterator( ++(*this) ); }
          //! Returns true if iterators are at same position.
          bool operator==(  const_iterator const& _b ) const 
            { return oiterator_ == _b.oiterator_; }
          //! Returns false if iterators are at same position.
          bool operator!=(  const_iterator const& _b ) const
            { return oiterator_ != _b.oiterator_; }
      
        private:
          //! initializes to first.
          const_iterator begin_( Bases<T_STRUCTURE> :: t_Structure const& _str );
          //! initializes to first.
          const_iterator end_( Bases<T_STRUCTURE> :: t_Structure const& _str );
          //! Value to return when dereferenced.
          Basis val_;
          //! origin iterator.
          typename Bases<T_STRUCTURE>::Origin oiterator_;
          //! origin iterator end.
          typename Bases<T_STRUCTURE>::Origin oiterator_end_;
          //! x iterator.
          typename Bases<T_STRUCTURE>::Xcoord xiterator_;
          //! x iterator end. 
          typename Bases<T_STRUCTURE>::Xcoord xiterator_end_;
          //! y iterator.
          typename Bases<T_STRUCTURE>::Ycoord yiterator_;
          //! y iterator.
          typename Bases<T_STRUCTURE>::Ycoord yiterator_end_;
          //! Size of the structure.
          size_t N_;
      };
      
  } // end of atomic_potential namespace.
} // namespace LaDa

#include "bases.impl.h"

#endif
