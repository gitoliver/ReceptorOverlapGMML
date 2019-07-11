#include "../includes/bead_residues.h"


AtomVector beads::Add_Beads_To_Protein(MolecularModeling::Assembly &assembly)
{
    AtomVector protein_beads;
    ResidueVector protein_residues = assembly.GetAllProteinResiduesOfAssembly();
    for (ResidueVector::iterator it1 = protein_residues.begin(); it1 != protein_residues.end(); ++it1)
    {
        Residue *residue = *it1;
        AtomVector atoms = residue->GetAtoms();
        for (AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); ++it2)
        {
            Atom *atom = *it2;
            // Okay so all CA atoms and then I've picked out others by hand in VMD that will cover the sidechains
            if (atom->GetName().compare("CA")==0)
            {
                //std::cout << "Adding bead to protein " << residue->GetId() << std::endl;
                Atom* bead_atom = new Atom(residue, "mfat", atom->GetCoordinates());
                residue->AddAtom(bead_atom);
                protein_beads.push_back(bead_atom);
            }
            else if ( (atom->GetName().compare("NZ")==0) ||
                      (atom->GetName().compare("CZ")==0) ||
                      (atom->GetName().compare("NE2")==0) ||
                      (atom->GetName().compare("OD1")==0) ||
                      (atom->GetName().compare("SD")==0)  ||
                      ( (atom->GetName().compare("CE2")==0) && residue->GetName().compare("TRP")==0 ) ||
                      ( (atom->GetName().compare("CD1")==0) && ( residue->GetName().compare("LEU")==0 || residue->GetName().compare("ILE")==0 ) ) ||
                      ( (atom->GetName().compare("CD")==0) && residue->GetName().compare("GLU")==0 )
                    )
            { // sfats should move when a chi1, chi2 is moved, so make sure they are connected to something for the SetDihedral function to move them.
                Atom* bead_atom = new Atom(residue, "sfat", atom->GetCoordinates());
                residue->AddAtom(bead_atom);
                protein_beads.push_back(bead_atom);
            }
        }
    }
    return protein_beads;
}

AtomVector beads::Add_Beads_To_Glycan(MolecularModeling::Assembly *assembly)
{
    AtomVector glycan_beads;
    ResidueVector glycan_residues = assembly->GetResidues();
    for (ResidueVector::iterator it2 = glycan_residues.begin(); it2 != glycan_residues.end(); ++it2)
    {
        Residue *residue = *it2;
        if ( residue->GetName().at(1) != 'S' ) // don't add a bead to a sialic acid in this section (see "else" below). Middle character of resname is always S for sialic acid.
       // if ( residue->GetName().compare(1, 2, "SA") != 0) // don't add one to a sialic acid (see below)
        {
            // std::cout << (*resi_iter)->GetName() << "\tG\t" << (*resi_iter)->CheckIfProtein() << endl;
           // std::cout << "Adding bead to self glycan " << residue->GetId() << std::endl;
            Atom* bead_atom = new Atom(residue, "gfat", residue->GetGeometricCenter());
            residue->AddAtom(bead_atom);
            glycan_beads.push_back(bead_atom);
            //Bond bead_atom to any other atom in residue so when glycan is moved, bead_atom moves too.
            Atom *any_atom = residue->GetAtoms().at(0); // 0 is arbitrary, any atom would do.
            //std::cout << "Blow here?" << any_atom->GetId() << std::endl;
            any_atom->GetNode()->AddNodeNeighbor(bead_atom);
            AtomVector temp = {any_atom};
            AtomNode *node = new AtomNode(); // DELETE IS FOR LOSERS.
            bead_atom->SetNode(node);
            bead_atom->GetNode()->SetNodeNeighbors(temp);
            AtomVector atoms = residue->GetAtoms();
            for (AtomVector::iterator it3 = atoms.begin(); it3 != atoms.end(); ++it3)
            {
                Atom *atom = *it3;
                if ( (atom->GetName().compare("C2N") == 0) || (atom->GetName().compare("C6") == 0) )
                {
                    bead_atom = new Atom(residue, "gfat", atom->GetCoordinates().at(0));
                    residue->AddAtom(bead_atom);
                    glycan_beads.push_back(bead_atom);
                    any_atom->GetNode()->AddNodeNeighbor(bead_atom);
                    AtomNode *node1 = new AtomNode(); // DELETE IS FOR LOSERS.
                    bead_atom->SetNode(node1);
                    bead_atom->GetNode()->SetNodeNeighbors(temp);
                }
            }
        }
        else // if it is sialic acid
        {
            AtomVector atoms = residue->GetAtoms();
            for (AtomVector::iterator it3 = atoms.begin(); it3 != atoms.end(); ++it3)
            {
                Atom *atom = *it3;
                if ( (atom->GetName().compare("C2") == 0) || (atom->GetName().compare("N5") == 0) || (atom->GetName().compare("C8") == 0) )
                {
                    Atom* bead_atom = new Atom(residue, "gfat", atom->GetCoordinates().at(0));
                    residue->AddAtom(bead_atom);
                    glycan_beads.push_back(bead_atom);
                    Atom *any_atom = residue->GetAtoms().at(0);
                    any_atom->GetNode()->AddNodeNeighbor(bead_atom);
                    AtomVector temp = {any_atom};
                    AtomNode *node = new AtomNode(); // DELETE IS FOR LOSERS.
                    bead_atom->SetNode(node);
                    bead_atom->GetNode()->SetNodeNeighbors(temp);
                }
            }
        }
    }
    return glycan_beads;
}

void beads::Remove_Beads(MolecularModeling::Assembly &assembly)
{
    ResidueVector all_residues = assembly.GetAllResiduesOfAssembly();
    for (ResidueVector::iterator it1 = all_residues.begin(); it1 != all_residues.end(); ++it1)
    {
        Residue *residue = *it1;
        AtomVector atoms = residue->GetAtoms();
        for (AtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); ++it2)
        {
            Atom *atom = *it2;
            if (atom->GetName().find("fat")==1)
            {
                residue->RemoveAtom(atom);
            }
        }
    }
}

double beads::Calculate_bead_overlaps(AtomVector &atomsA, AtomVector &atomsB)
{
    double radius = 3.0; //Using same radius for all beads.
    double distance = 0.0, overlap = 0.0;
  //  std::cout << "About to check " << atomsA.size() << " atoms vs " << atomsB.size() << " atoms" << std::endl;
    for(AtomVector::iterator it1 = atomsA.begin(); it1 != atomsA.end(); ++it1)
    {
        Atom *atomA = *it1;
        for(AtomVector::iterator it2 = atomsB.begin(); it2 != atomsB.end(); ++it2)
        {
            Atom *atomB = *it2;
            if ( (atomA->GetCoordinates().at(0)->GetX() - atomB->GetCoordinates().at(0)->GetX()) < (radius * 2) ) // This is faster than calulating distance, and rules out tons of atom pairs.
            {
                distance = atomA->GetDistanceToAtom(atomB);
                if ( ( distance < (radius + radius) ) && ( distance > 0.0 ) ) //Close enough to overlap, but not the same atom
                {
                    overlap += gmml::CalculateAtomicOverlaps(atomA, atomB, radius, radius); // This calls the version with radius values
                }
            }
        }
    }
    //return (overlap );
    return (overlap / gmml::CARBON_SURFACE_AREA); //Normalise to area of a buried carbon
}





