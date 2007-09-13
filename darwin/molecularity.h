//
//  Version: $Id$
//
#ifndef _MOLECULARITY_H_
#define _MOLECULARITY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <vff/functional.h>
#include <pescan_interface/interface.h>
#include <lamarck/structure.h>
#include <opt/opt_function_base.h>
#include <opt/opt_minimize_gsl.h>
#include <opt/types.h>

#include "pescan.h"
#include "evaluator.h"
#include "concentration.h"
#include "functors.h"
#include "individual.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace Molecularity
{
  // Exact same object as BandGap, except for the two friends
  struct Object : public BandGap::Object
  {
    friend void operator<<(Ising_CE::Structure &_str, const Object &_o)
    friend operator<<(Object &_o, const Ising_CE::Structure &_c)
    Object() : BandGap::Object() {}
    Object(const BandGap::Object &_c) : BandGap::Object(_c) {}
    Object(const Object &_c) : BandGap::Object(_c) {}
    ~Object() {};
  };


  class Evaluator : public BandGap::Evaluator< Individual::Types<BandGap::Object>::Scalar >
  {
    protected:
      class BandGap::Evaluator< Individual::Types<BandGap::Object>::Scalar > t_Base;
    public:
      typedef  t_Base :: t_Individual t_Individual;

    public:
      Evaluator() : t_Base() {}
      ~Evaluator() {}


      bool initialize( t_Individual &_indiv );
  }

} // namespace BandGap


#include "molecularity.impl.h"

#endif // _MOLECULARITY_H_
