#ifndef _FUNCTIONAL_BUILDER_H_
#define _FUNCTIONAL_BUILDER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <stdio.h>

#include <opt/opt_functors.h>
#include <tinyxml/tinyxml.h>

#include "atat/vectmac.h"

#include "polynome.h"
#include "cluster.h"
#include "structure.h"
#include "constituent_strain.h"
#include "lattice.h"

#ifdef _MPI 
  #include "mpi/mpi_object.h"
#endif

namespace VA_CE {

  using Ising_CE::Constituent_Strain;
  using Ising_CE::Atat_Structure;

  class Functional_Builder
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<Functional_Builder> ( Functional_Builder& );
#endif
    public:
      typedef function::Plus< Polynome, Constituent_Strain > t_VA_Functional;

    protected:
      static Ising_CE::Lattice*              lattice;
      static std::vector<Ising_CE::Cluster>* clusters;
      static Constituent_Strain*             harmonics;
      
    // constructor, destructor
    public:
      Functional_Builder(){};
      virtual ~Functional_Builder();

      virtual bool Load( const TiXmlElement &handle );

    // generators
    public: 
      void add_equivalent_clusters();
      bool generate_functional(const Ising_CE::Structure &str, t_VA_Functional * const functional);

  };

  
} // end of namespace VA_CE
#endif
