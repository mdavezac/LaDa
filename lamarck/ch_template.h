//
//  Version: $Id$
//
#ifndef _CH_TEMPLATE_H_
#define _CH_TEMPLATE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include <tinyxml/tinyxml.h>

#include "opt/types.h"
#include "atat/vectmac.h"

#include "structure.h"

namespace VA_CE {

  class CH_Template
  {
    public:
      CH_Template() {};
      virtual ~CH_Template() {};
      virtual bool add_structure( const types::t_real _energy,
                                  const Ising_CE::Structure &_structure ) = 0;

      // required lada::Fitness_Function behaviors
      virtual types::t_real evaluate() const = 0;
      virtual types::t_real evaluate(const types::t_real _x) = 0;
      virtual types::t_real evaluate_gradient() = 0;
      virtual types::t_real evaluate_gradient(const types::t_real _x) = 0;

      // other
      virtual bool Load(const TiXmlElement &_element) = 0;
      virtual void print_out (std::ostream &stream, types::t_int print_what) const = 0;
      virtual void print_xml( TiXmlElement* const node ) const = 0;
  };


} // namespace VA_CE
#endif
