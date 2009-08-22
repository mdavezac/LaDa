//
//  Version: $Id$
//
#ifndef LADA_ATOMIC_POTENTIAL_COLLAPSE_FITTING_SETS_H_
#define LADA_ATOMIC_POTENTIAL_COLLAPSE_FITTING_SETS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <utility>

#include <opt/debug.h>

#include "../numeric_types.h"


namespace LaDa
{
// //! \cond
// namespace Crystal
// {
//   template<class T_TYPE> class TStructure;
// }
// //! \endcond

  namespace atomic_potential
  {
    // Forward declaration, 
    //! \cond
    class Representation;
    //! \endcond

    namespace collapse
    {
      //! \brief   Fitting structure, with coordinates, weights, and energies.
      //! \details Input sets for fitting are described by a number of
      //!          structures with associated weights and energies. 
      //!          In the case of atomic potential, a single structure is
      //!          described by a number of symmetric representation. This
      //!          class helps deal with the plurality of representations. 
      //!          Containers are arranged such that coordinates of the
      //!          representation are the outermost index. This unnatural
      //!          storage helps perform the alternating least-square fit.
      class FittingSet
      {
        protected:
          //! Type of the coordinates for a single representation.
          typedef std::vector<specie_type> t_AtomType;
          //! Type of container over coordinates for each  structure and their representations.
          typedef std::vector // container over coodinates
                  <
                    std::vector // for each structure in input set,
                    <
                      std::vector  // for each representation of each structure.
                      <
                        t_AtomType
                      >
                    >
                  > t_Coordinates;
          //! Type of the containers for weights.
          typedef std::vector // over structures.
                  <
                    std::pair
                    <
                      numeric_type, // fitting weight.
                      std::vector<numeric_type> // weight of each representation.
                    >
                  > t_Weights;
          //! Energies of each structure.
          typedef std::vector<numeric_type> t_Energies;

        public:
          //! An iterator over the structures.
          class str_iterator;
          //! Constructor.
          FittingSet() {}
          //! Copy Constructor.
          FittingSet    (FittingSet const &_c)
                      : coordinates_(_c.coordinates_),
                             energies_(_c.energies_), 
                             weights_(_c.weights_){}
          //! Changes the weight of nth structure added to representations.
          numeric_type change_weight( size_t _i, numeric_type _w)
            { numeric_type const old(weights_[_i].first); weights_[_i].first=_w; return old; } 
          //! Adds a structure to the representations.
          void add(Representation const &_representation,
                   numeric_type _energy, numeric_type _weight);
          //! returns iterator to start of structures, for variable \a _i.
          str_iterator begin( size_t _i ) const;
          //! returns iterator to end of structures, for variable \a _i.
          str_iterator end( size_t _i ) const;

        protected:
          //! Coordinates for each structure and representation.
          t_Coordinates coordinates_;
          //! Weight of each structure and represenation.
          t_Energies energies_;
          //! Weight of the represenations.
          t_Weights weights_;
      };

      //! True if iterators are at same position.
      bool operator==( FittingSet::str_iterator const& _a,
                       FittingSet::str_iterator const& _b );

      //! \brief Helps iterate over representations.
      //! \warning This special iterator cannot be dereferenced directly. 
      class FittingSet::str_iterator
      {
        friend bool operator==( str_iterator const& _a, str_iterator const& _b );
        friend class FittingSet;

        public:
          //! Iterator over representations.
          class rep_iterator;
          //! Constructor.
          str_iterator() {}
          //! Copy Constructor.
          str_iterator   (str_iterator const &_c) 
                       : i_(_c.i_),
                         i_coordinates_(_c.i_coordinates_),
                         i_energy_(_c.i_energy_),
                         i_weight_(_c.i_weight_) {}

          //! Increments iterator. Returns true if not end of container.
          void operator++() {  ++i_coordinates_; ++i_weight_; ++i_energy_; }

          //! Return current weight.
          numeric_type weight() const { return i_weight_->first; }
          //! Return current energy.
          numeric_type energy() const { return *i_energy_; }
          //! Returns iterator over representations.
          rep_iterator begin() const;
          //! Returns iterator over representations.
          rep_iterator end() const;

        protected:
          //! Current variable.
          size_t i_;
          //! Iterator over input structures for specific variable.
          t_Coordinates::value_type::const_iterator i_coordinates_;
          //! Iterator over energies.
          t_Energies::const_iterator i_energy_;
          //! Iterator over fitting weights of structures for specific variable.
          t_Weights::const_iterator i_weight_;
      };

      //  True if iterators are at same position.
      inline bool operator==( FittingSet::str_iterator const& _a,
                              FittingSet::str_iterator const& _b )
      {
        LADA_ASSERT( _a.i_ == _b.i_, "Inequivalent iterators.\n");
        return _a.i_energy_ != _b.i_energy_; 
      }
      //! False if iterators are at same position.
      inline bool operator!=( FittingSet::str_iterator const& _a,
                              FittingSet::str_iterator const& _b )
        { return not (_a == _b); }

      //! True if iterators are at same position.
      bool operator==( FittingSet::str_iterator::rep_iterator const& _a,
                       FittingSet::str_iterator::rep_iterator const& _b );
      
      class FittingSet::str_iterator::rep_iterator
      {
        friend bool operator==( rep_iterator const& _a, rep_iterator const& _b );
        friend class FittingSet::str_iterator;
        public:
          typedef t_AtomType :: const_iterator coordinate_iterator;
          //! Constructor.
          rep_iterator() {}
          //! Copy Constructor.
          rep_iterator   (rep_iterator const &_c) 
                       : i_coordinates_(_c.i_coordinates_),
                         i_weight_(_c.i_weight_)  {}
          //! Returns iterator to coordinates.
          coordinate_iterator begin() const { return i_coordinates_->begin(); }
          //! Returns iterator to coordinates.
          coordinate_iterator end() const { return i_coordinates_->end(); }
          //! Returns the weight.
          t_Weights::value_type::second_type::value_type weight() const { return *i_weight_; }
          //! Increments iterators.
          void operator++() { ++i_coordinates_; ++i_weight_; }
          
        protected:
          //! Iterator over input structures for specific variable.
          t_Coordinates::value_type::value_type::const_iterator i_coordinates_;
          //! Iterator over fitting weights of structures for specific variable.
          t_Weights::value_type::second_type::const_iterator i_weight_;
      };

      //  True if iterators are at same position.
      inline bool operator==( FittingSet::str_iterator::rep_iterator const& _a,
                              FittingSet::str_iterator::rep_iterator const& _b )
        { return _a.i_coordinates_ != _b.i_coordinates_; }
      //! False if iterators are at same position.
      inline bool operator!=( FittingSet::str_iterator::rep_iterator const& _a,
                              FittingSet::str_iterator::rep_iterator const& _b )
        { return not (_a == _b); }
    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
#endif

