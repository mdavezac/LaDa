#include "LaDaConfig.h"
#include<iostream>
#include<string>

#include <opt/debug.h>
#include "../symmetry_operator.h"

using namespace std;
bool check_matrix(LaDa::math::SymmetryOperator const &_sym, theta)
int main()
{
  using namespace LaDa;
  using namespace LaDa::math;
  
  SymmetryOperator symop();
  types::t_real theta = 0.8;
  symop.set_rotation(1, 0, 0)
                    (0, std::cos(theta), -std::sin(theta))
                    (0, std::sin(theta),  std::cos(theta));
  LADA_DOASSERT(symop.is_symmetry(), "Failed on symmetry.\n")

  namespace lns = LaDa::load_n_save;
  Atom< LADA_TYPE > atom;
  atom.pos = math::rVector3d(0.5,-0.5,0.5);
  atom.freeze = Atom< LADA_TYPE >::frozen::X;
  LADA_INIT_TYPE;

  // print to xml
  lns::save::Save saver;
  boost::shared_ptr< lns::tree::Base > other(saver(lns::ext(atom)));
  std::ostringstream sstr;
  lns::xml::print(sstr, *other);
  std::string xmlstring = sstr.str();
  boost::algorithm::erase_all(xmlstring, "\n");
  boost::algorithm::trim(xmlstring);
  LADA_DOASSERT(xmlstring == LADA_XML, "Error in output.");

  atom.pos = math::rVector3d(0,0,0);
  atom.freeze = Atom< LADA_TYPE >::frozen::NONE;
  LADA_CLEAR_TYPE;
  other = lns::xml::parse( sstr.str() );
  lns::load::Load loader;
  bool result = loader( *other, lns::ext(atom) );
  LADA_DOASSERT((atom.pos - math::rVector3d(0.5,-0.5,0.5)).squaredNorm() < 1e-12, "Could not reload pos.\n");
  LADA_DOASSERT(atom.freeze == Atom< LADA_TYPE >::frozen::X, "Could not reload freeze.\n");
  LADA_DOASSERT_TYPE

  return 0;
}
