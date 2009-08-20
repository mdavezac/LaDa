//
//  Version: $Id$
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
  namespace atomic_potential
  {
    struct Basis 
    {
      atat::rVector3d origin;
      atat::rVector3d x;
      atat::rVector3d y;
      atat::rVector3d z;
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

    //! Comparison of origin iterators.
    template<class T_STRUCTURE>
      bool operator==( Bases<T_STRUCTURE>::Origin const&, Bases<T_STRUCTURE>::Origin const& );

    template<class T_STRUCTURE>
      class Bases::Origin
      {
        friend class Xcoord;
        friend class Ycoord;
        friend bool operator==( Origin const&, Origin const& );
        public:
          //! Return on deref.
          typedef atat::rVector3d const& value_type;
          //! Return on deref.
          typedef atat::rVector3d const* pointer_type;
          //! Constructor.
          Origin() {}
          //! Copy Constructor.
          Origin( Origin const& _c ) : structure_(_c.structure_), iterator_(_c.iterator_) {};
  
          //! Deref.
          value_type operator*() const { return iterator_->pos; }
          //! Deref.
          pointer_type operator->() const { return &(iterator_->pos); }
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
    
    template<class T_STRUCTURE>
      inline bool operator==( Bases<T_STRUCTURE>::Origin const& _a,
                              Bases<T_STRUCTURE>::Origin const& _b )
        { return _a.iterator_ == _b.iterator_; }
    template<class T_STRUCTURE>
      inline bool operator!=( Bases<T_STRUCTURE>::Origin const& _a,
                              Bases<T_STRUCTURE>::Origin const& _b )
      { return not(_a == _b); }

    //! Comparison of x iterators.
    template<class T_STRUCTURE>
      bool operator==( Bases<T_STRUCTURE>::Xcoord const&, Bases<T_STRUCTURE>::Xcoord const& );
    
    
    //! Type of the iterator over x axis.
    class Bases<T_STRUCTURE> :: Xcoord 
    {
      friend class Ycoord;
      friend bool operator==( Xcoord const&, Xcoord const& );
      //! Type of the vector of neighbors.
      typedef std::vector< atat::rVector3d > t_Neighs;
      public:
        typedef atat::rVector3d const& value_type;
        typedef atat::rVector3d const* pointer_type;
    
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

    inline bool operator==( Bases<T_STRUCTURE>::Xcoord const& _a,
                            Bases<T_STRUCTURE>::Xcoord const& _b )
      { return _a.iterator_ == _b.iterator_; }
    inline bool operator!=( Bases<T_STRUCTURE>::Xcoord const& _a,
                            Bases<T_STRUCTURE>::Xcoord const& _b )
      { return not(_a == _b); }

    //! Comparison of y iterators.
    bool operator==( Bases<T_STRUCTURE>::Ycoord const& _a,
                     Bases<T_STRUCTURE>::Ycoord const& _b );

    class Bases<T_STRUCTURE> :: Ycoord
    {
      friend bool operator==( Ycoord const&, Ycoord const& );
      public:
        typedef atat::rVector3d const& value_type;
        typedef atat::rVector3d const* pointer_type;

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
        size_t size() const { return equivs_ ? equivs_->size(): 0; }

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

    inline bool operator==( Bases<T_STRUCTURE>::Ycoord const& _a,
                            Bases<T_STRUCTURE>::Ycoord const& _b )
      { return _a.iterator_ == _b.iterator_; }
    inline bool operator!=( Bases<T_STRUCTURE>::Ycoord const& _a,
                            Bases<T_STRUCTURE>::Ycoord const& _b )
      { return not(_a == _b); }

    //! Comparison of basis iterators.
    bool operator==( Bases<T_STRUCTURE>::const_iterator const& _a,
                     Bases<T_STRUCTURE>::const_iterator const& _b );

    class Bases<T_STRUCTURE>::const_iterator
    {
      friend class Bases<T_STRUCTURE>;
      friend bool operator==( const_iterator const&, const_iterator const& );
      public:
        typedef Basis const& value_type;
        typedef Basis const* pointer_type;

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
        value_type operator*() const { return val_; }
        //! Deref.
        pointer_type operator->() const { return &val_; }
        //! Pre-increment operator.
        const_iterator &operator++();
        //! Post-increment operator.
        const_iterator operator++(int) { return const_iterator( ++(*this) ); }

      private:
        //! initializes to first.
        const_iterator begin_( Bases<T_STRUCTURE> :: t_Structure const& _str );
        //! initializes to first.
        const_iterator end_( Bases<T_STRUCTURE> :: t_Structure const& _str );
        //! Value to return when dereferenced.
        Basis val_;
        //! origin iterator.
        Bases<T_STRUCTURE>::Origin oiterator_;
        //! origin iterator end.
        Bases<T_STRUCTURE>::Origin oiterator_end_;
        //! x iterator.
        Bases<T_STRUCTURE>::Xcoord xiterator_;
        //! x iterator end. 
        Bases<T_STRUCTURE>::Xcoord xiterator_end_;
        //! y iterator.
        Bases<T_STRUCTURE>::Ycoord yiterator_;
        //! y iterator.
        Bases<T_STRUCTURE>::Ycoord yiterator_end_;
        //! Size of the structure.
        size_t N_;
    };

    inline bool operator==( Bases<T_STRUCTURE>::const_iterator const& _a,
                            Bases<T_STRUCTURE>::const_iterator const& _b )
      { return _a.oiterator_ == _b.oiterator_; }
    inline bool operator!=( Bases<T_STRUCTURE>::const_iterator const& _a,
                            Bases<T_STRUCTURE>::const_iterator const& _b )
      { return not(_a == _b); }
  } // end of atomic_potential namespace.
} // namespace LaDa

#endif
