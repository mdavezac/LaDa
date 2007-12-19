//
//  Version: $Id$
//
#ifndef _FUNCTIONAL_BUILDER_H_
#define _FUNCTIONAL_BUILDER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <stdio.h>

#include <opt/function_functors.h>
#include <tinyxml/tinyxml.h>

#include <atat/vectmac.h>

#include "polynome.h"
#include "cluster.h"
#include "structure.h"
#include "constituent_strain.h"
#include "lattice.h"

#ifdef _MPI 
  #include <mpi/mpi_object.h>
#endif

//! \brief Nonsense namespace which should become %CE (probably)
namespace VA_CE {

  using Ising_CE::Constituent_Strain;
  using Ising_CE::Atat_Structure;

  //! \brief This creates a cell-shape specialized functional out of a set of clusters.
  //! \details The cluster expansion formalism is meant to be true for the
  //!          full infinite lattice. Unfortunately, within such a paradigm,
  //!          applying the cluster expansion to a particular structure is not
  //!          effective. This class transforms a set of cluster interactions
  //!          and harmonics to a cell-shape specialized functional constituted
  //!          of a polynomial and a constituent strain instance.
  class Functional_Builder
  {
#ifdef _MPI
    //! \cond
    friend bool mpi::BroadCast::serialize<Functional_Builder> ( Functional_Builder& );
    //! \endcond
#endif
    public:
      //! Type of the specialized functional
      typedef function::Plus< Polynome, Constituent_Strain > t_VA_Functional;

    protected:
      //! Pointer to the lattice 
      static Ising_CE::Lattice*              lattice;
      //! Pointer to the set of clusters
      static std::vector<Ising_CE::Cluster>* clusters;
      //! Pointer to the harmonics
      static Constituent_Strain*             harmonics;
      
    // constructor, destructor
    public:
      //! Constructor
      Functional_Builder() {}
      //! Copy Constructor
      Functional_Builder   (const Functional_Builder &_c)  {}
      //! Destructor
      virtual ~Functional_Builder();

      //! Loads a %Cluster Expansion functional from XML
      virtual bool Load( const TiXmlElement &handle );

    public: 
      //! \brief Finds all symmetrically equivalent clusters of the clusters in
      //!        Functional_Builder::clusters.
      //! \details In effect, if those clusters in Functional_Builder::Clusters
      //!          a generating family prior to call, then after call,
      //!          Functional_Builder::Clusters forms a group with respect to
      //!          every potin-symmetry operation of the lattice,
      void add_equivalent_clusters();
      //! \brief Creates a specialized cluster expansion functional for structure \a str.
      bool generate_functional(const Ising_CE::Structure &str, t_VA_Functional * const functional);

  };

  
} // end of namespace VA_CE
#endif
