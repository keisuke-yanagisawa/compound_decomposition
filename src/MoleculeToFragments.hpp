#include "common.hpp"
#include "Fragment.hpp"
#include "Molecule.hpp"
#ifndef MOLECULE_TO_FRAGMENTS_H_
#define MOLECULE_TO_FRAGMENTS_H_

namespace fragdock {
  std::vector<Fragment> getFragments(const Molecule &mol, 
                                     int max_ring_size=-1,
                                     bool merge_solitary=true);
}
#endif
