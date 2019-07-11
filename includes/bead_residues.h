#ifndef BEAD_RESIDUES_H
#define BEAD_RESIDUES_H

#include "/home/oliver/Programs/gems/gmml/includes/gmml.hpp" // can just do gmml.hpp, but for qt autocomplete path is better
//#include "selections.h"

using namespace MolecularModeling;

namespace beads
{

void Remove_Beads(MolecularModeling::Assembly &assembly);
AtomVector Add_Beads_To_Glycan(MolecularModeling::Assembly *assembly);
AtomVector Add_Beads_To_Protein(MolecularModeling::Assembly &assembly);
double Calculate_bead_overlaps(AtomVector &atomsA, AtomVector &atomsB);

}
#endif // BEAD_RESIDUES_H

