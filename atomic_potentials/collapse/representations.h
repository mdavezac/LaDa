//
//  Version: $Id$
//
#ifndef LADA_ATOMIC_POTENTIAL_COLLAPSE_WEIGHTS_H_
#define LADA_ATOMIC_POTENTIAL_COLLAPSE_WEIGHTS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa
{
  namespace atomic_potential
  {
    namespace collapse
    {
      //! \brief   Weight of a fitting structure.
      //! \details Input sets for fitting are described by a number of
      //!          structures with associated weights and energies. 
      //!          In the case of atomic potential, a single structure is
      //!          described by a number of symmetric representation. This
      //!          class helps deal with the plurality of representations. 
      //!          Containers are arranged such that coordinates of the
      //!          representation are the outermost index. This unnatural
      //!          storage helps perform the alternating least-square fit.
      class Representations
      {
        protected:
          //! Type of container over coordinates for each  structure and their representations.
          typedef std::vector // container over coodinates
                  <
                    std::vector // for each structure in input set,
                    <
                      std::vector  // for each representation of each structure.
                      <
                        types::t_real
                      >
                    >
                  > t_Coordinates;
          //! Type of the containers for weights.
          typedef std::vector // over structures.
                  <
                    std::pair
                    <
                      types::t_real, // fitting weight.
                      std::vector<types::t_real> // weight of each representation.
                    >
                  > t_Weights;
          //! Energies of each structure.
          typedef std::vector<types::t_real> t_Energies;

        public:
          //! An iterator over the representations.
          class const_iterator;
          //! Constructor.
          Representations() {}
          //! Copy Constructor.
          Representations    (Representations const &_const)
                           : coordinates_(_c.coordinates_),
                             energies_(_c.energies_), 
                             weights_(_c.weights_){}
          //! Changes the weight of nth structure added to representations.
          void change_weight( size_t _i, types::t_real _w)
            { weights_[_i].first = _w; } 
          //! Adds a structure to the representations.
          void add( Crystal::Structure &_structure, types::t_real _w );

        protected:
          //! Coordinates for each structure and representation.
          t_Coordinates coordinates_;
          //! Weight of each structure and represenation.
          t_Energies energies_;
          //! Weight of the represenations.
          t_Weights weights_;
      };

      //! \brief Helps iterate over representations.
      //! \warning This special iterator cannot be dereferenced directly. 
      class Representations::const_iterator
      {
        public:
          //! Constructor.
          const_iterator() {}
          //! Copy Constructor.
          const_iterator   (const_iterator const &_c) 
                         : i_str_coors_(_c.i_str_coors_),
                           i_rep_coors_(_c.i_rep_coors_),
                           i_rep_coors_end_(_c.i_rep_coors_end_),
                           i_str_weight_(_c.i_str_weights_),
                           i_rep_weight_(_c.i_rep_weight_)
#                          ifdef LADA_DEBUG
                             , i_rep_weight_end_(_c.i_rep_weight_end_)
#                          endif
                           {}

          const_iterator& operator++() 
          {
            if( i_str_coors_ == i_str_coors_end_ ) return *this;
            ++i_rep_coors;
            if( i_rep_coors_ == i_rep_coors_end_)
            {
              LADA_ASSERT( i_rep_weight_ == i_rep_weight_end_, "Incoherent container.\n" );
              ++i_str_coors_;
            }

          }


        protected:
          //! Iterator over input structures for specific variable.
          t_Coordinates::value_type::const_iterator i_str_coors_;
          //! Iterator over input structures for specific variable.
          t_Coordinates::value_type::const_iterator i_str_coors_end_;
          //! Iterator over representation structures for specific variable.
          t_Coordinates::value_type::value_type::const_iterator i_rep_coors_;
          //! Iterator over representation structures for specific variable.
          t_Coordinates::value_type::value_type::const_iterator i_rep_coors_end_;
          //! Iterator over energies.
          t_Energies::const_iterator i_energs;
          //! Iterator over fitting weights of structures for specific variable.
          t_Weights::const_iterator i_str_weight_;
          //! Iterator over representation weights of structures for specific variables.
          t_Weights::const_iterator i_rep_weight_;
#         ifdef LADA_DEBUG
            //! Iterator over representation weights of structures for specific variables.
            t_Weights::const_iterator i_rep_weight_end_;
#         endif
      };
    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
#endif

