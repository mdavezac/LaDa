#ifndef _FUNCTIONAL_BUILDER_H_
#define _FUNCTIONAL_BUILDER_H_

#include "LaDaConfig.h"

#include <string>
#include <stdio.h>

#include <opt/function_functors.h>
#include <tinyxml/tinyxml.h>
#include <mpi/mpi_object.h>
#include <crystal/structure.h>
#include <crystal/lattice.h>

#include "polynome.h"
#include "cluster.h"
#include "constituent_strain.h"


namespace LaDa
{

  //! Holds %Cluster Expansion related stuff.
  namespace CE 
  {

    //! \brief This creates a cell-shape specialized functional out of a set of clusters.
    //! \details The cluster expansion formalism is meant to be true for the
    //!          full infinite lattice. Unfortunately, within such a paradigm,
    //!          applying the cluster expansion to a particular structure is not
    //!          effective. This class transforms a set of cluster interactions
    //!          and harmonics to a cell-shape specialized functional constituted
    //!          of a polynomial and a constituent strain instance.
    template< class T_HARMONIC >
    class Builder
    {
      public:
        //! Type of the harmonic used.
        typedef T_HARMONIC t_Harmonic;
        //! Type of the constituent strain.
        typedef ConstituentStrain :: Functional< t_Harmonic > t_CS;
        //! Type of the chemical functional.
        typedef Polynome t_Chemical;
        //! Type of the specialized functional
        typedef function::Plus< t_Chemical, t_CS > t_VA_Functional;
        //! Type of the clusters.
        typedef Cluster t_Cluster;
        //! Type of the vector of clusters.
        typedef std::vector< t_Cluster > t_Clusters;

      protected:
        //! Pointer to the lattice 
        static Crystal::Lattice* lattice;
        //! Pointer to the set of clusters
        t_Clusters*        clusters;
        //! Pointer to the harmonics
        static t_CS*              harmonics;
        
      // constructor, destructor
      public:
        //! Constructor
        Builder() {}
        //! Copy Constructor
        Builder   (const Builder &_c)  {}
        //! Destructor
        virtual ~Builder();

        //! Loads a %Cluster Expansion functional from XML
        virtual bool Load( const TiXmlElement &handle );

      public: 
        //! \brief Finds all symmetrically equivalent clusters of the clusters in
        //!        Builder::clusters.
        //! \details In effect, if those clusters in Builder::Clusters
        //!          a generating family prior to call, then after call,
        //!          Builder::Clusters forms a group with respect to
        //!          every symmetry operation of the lattice,
        void add_equivalent_clusters() { CE::add_equivalent_clusters(*lattice, *clusters); }
        //! \brief Creates a specialized cluster expansion functional for structure \a str.
        bool generate_functional(const Crystal::Structure &str,
                                 t_VA_Functional * const functional) const;
        //! \brief Builds and returns a pair of chemical and strain functionals.
        //! \details The functional are specialized for the structure on input.
        //! \warning the variables of the functionals are not set to anything.
        //!          see function::Base for details.
        std::pair<t_Chemical*, t_CS*> generate_functional(const Crystal::Structure &_str ) const;
    };

    
  } // end of namespace CE

} // namespace LaDa

#include "functional_builder.impl.h"

#endif
