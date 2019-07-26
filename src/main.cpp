#include <iostream>
#include <fstream>
#include <iostream>
#include <cstdlib> // for exit()
// Includes for directory reading
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>     /* getenv */

#include "/home/oliver/Programs/gems/gmml/includes/gmml.hpp" // can just do gmml.hpp, but for qt autocomplete path is better
//#include "gmml.hpp"
#include "../includes/io.h"
#include "../includes/residue_name_conversion_map.h"
#include "../includes/bead_residues.h"

//typedef std::vector<GeometryTopology::Coordinate> CoordinateVector;
double GetAngleBetweenVectors(Eigen::Vector3d v1, Eigen::Vector3d v2, bool returnDegrees = false);
Atom* FindCAAtomOfResidueWithNumber(ResidueVector residues, std::string query_residue_number);
AssemblyVector loadPDBFilesFromDirectory(std::string directory);
AssemblyVector loadPDBFilesFromDirectory(std::string directory, MolecularModeling::Assembly combined_assembly);
MolecularModeling::Residue* seek_alycone_residue_connected_to_residue(MolecularModeling::Residue *residue, MolecularModeling::Assembly *assembly);
MolecularModeling::Residue* recursively_seek_alycone_residue(MolecularModeling::Atom *start_atom, bool *found);
double calculateInterStrepDistance(GeometryTopology::Coordinate *streptavidinCenterPoint, GeometryTopology::Coordinate *influenzaHAPoint1, GeometryTopology::Coordinate *influenzaHAPoint2);
//GeometryTopology::Coordinate getGeometricCenterOfCoordinates(CoordinateVector coordinates);

int main(int argc, char *argv[])
{

    /* Ok this is frankencode now. Originally adapted from Superimposing ERManI code, I wanted to make it general, but failed. The ERManI code
     * involved superimposing ERManI to a glycan and checking overlap. I needed to report the protein residue the glycan was attached to. So I had to
     * first find it and then convert it to HXB2 (HIV env) numbering. All that stuff I don't need but was attempting to remain flexible.
     * For this code, I have simulations (multiple snapshots/coordinate sets) of Strepavidin, which has four binding sites and the same biotinylated glycan
     * in each site. I want to do three things. First superimpose influenza HA to each glycan in each site and check if there overlaps between strep
     * and the HA (preventing binding). Second I want to measure the distance between the C1 of the Gal residue of the biotinylated glycan to the middle
     * of the strep molecule, as my collaborater is interested in the spacing required for high affinity binding (he knows this) and how high affinity
     * binding is being acheived (we can deduce from the spacing). Third, but related to the other two, I want to know the angle between the HA and the strep
     * as the HA is embedded in a viral coat in the actual experiment, and so if I predict that it can bind or am measureing distances from glycans, but
     * that shape wouldn't allow an embedded HA to bind then it shouldn't be counted. The first two questions could be solved with separate code. The third
     * question is really a pre-screen that I could do while doing the overlap calculation. Anyway I do them all in here for whatever reason that is
     * unclear to me now.


/* What do we want to know?
 * The overlap between ERManI and a glycoprotein after superimposition of ERManI to each of it's glycans, but not including the glycan we superimposed onto.
 * Also, want to know how the overlap changes over multiple snapshots
 * Also, want to differentiate between protein and glycan overlap.
 * Also, be able to report using HXB2 numbering for each site. Is a trimer so each site is repeated in the structure with unique residue numbers.
 *
 *
 * Load in pdb file(s) / trajectory of a glycosylated protein (Fuck it, just use an assembly vector).
 * load in pdb file of ERManI
 * load in pdb file of ligand of ERManI
 * Select out VMB only residues?
 * Superimpose ERManI to each VMB defined by residue number in inputs
 * Protein overlaps are easy, can select for protein only residues, but need overlaps with other glycans, and not current residue's glycan.
 *
 *
 * Reporting:
 * load in residue conversion map?
 * load in glycofrom conversion map?
 */

    if ( argc != 3 )
    {
        std::cout << "Usage: " << argv[0] << " relativePathToInputFolder/" << " input.txt\n";
        std::cout << "Example: " << argv[0] << " inputs/" << " input.txt\n";
        std::exit(1); // Take your ball and go home.
    }

    // Good comment.
    // Hardcoding as no time and this is too particular to be generalisable anyway. These first two points define the line that runs through the center
    // of the Strepavidin molecule.

    GeometryTopology::Coordinate streptavidinDisectingLinepoint2(40.1895, 21.084, 7.041);
    GeometryTopology::Coordinate streptavidinDisectingLinepoint1(6.096, 21.8635, -6.911);
    GeometryTopology::Coordinate streptavidinCenterPoint(23.143, 21.474, 0.065);

    GeometryTopology::Coordinate jola(1.0, -15.0, 1.0); // This is just me checking that this next function works

    // Now for the angle between strep and HA, need a vector for the step line and can calculate now as it is fixed.
    Eigen::Vector3d strepVector3D(streptavidinDisectingLinepoint1.GetX() - streptavidinDisectingLinepoint2.GetX(), streptavidinDisectingLinepoint1.GetY() - streptavidinDisectingLinepoint2.GetY(),streptavidinDisectingLinepoint1.GetZ() - streptavidinDisectingLinepoint2.GetZ() );
    Eigen::Vector3d strepVector3DReversed(streptavidinDisectingLinepoint2.GetX() - streptavidinDisectingLinepoint1.GetX(), streptavidinDisectingLinepoint2.GetY() - streptavidinDisectingLinepoint1.GetY(),streptavidinDisectingLinepoint2.GetZ() - streptavidinDisectingLinepoint1.GetZ() );

  //  std::cout << "Distance from jola to line is " << GeometryTopology::calculateDistanceFromPointToLineBetweenTwoPoints(jola, streptavidinDisectingLinepoint1, streptavidinDisectingLinepoint2) << "\n";

    typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;

    std::string workingDirectory = Find_Program_Working_Directory();
    std::string inputsDirectory = argv[1];
    std::string inputFile = argv[2];
    //************************************************//
    // Reading input file                             //
    //************************************************//

    StringVector target_residue_list, receptor_vector_residue_list;
    std::string buffer, snapshot_directory, template_pdb, ligand_pdb, receptor_pdb, residue_conversion_map;
    bool isResidueConversionMap = false, isReceptorVectorResidueList = false;
    std::ifstream infile (workingDirectory + "/" + inputsDirectory + "/" + inputFile);
    if (!infile)
    {
        std::cerr << "Uh oh, input file could not be opened for reading!" << std::endl;
        std::exit(1);
    }
    while (infile) // While there's still stuff left to read
    {
        std::string strInput;
        getline(infile, strInput);
        if(strInput == "Snapshot directory:")
            getline(infile, snapshot_directory);
        if(strInput == "Ligand PDB:")
            getline(infile, ligand_pdb);
        if(strInput == "Receptor PDB:")
            getline(infile, receptor_pdb);
        if(strInput == "Template PDB:")
            getline(infile, template_pdb);
        if(strInput == "Receptor Vector Residue Numbers:")
        {
            isReceptorVectorResidueList = true;
            getline(infile, buffer);
            while(buffer != "END")
            {
                receptor_vector_residue_list.push_back(buffer);
                getline(infile, buffer);
            }
            buffer.clear();
        }
        if(strInput == "Residue Conversion Map:")
        {
            isResidueConversionMap = true;
            getline(infile, residue_conversion_map);
        }
        if(strInput == "Target Residues:")
        {
            getline(infile, buffer);
            while(buffer != "END")
            {
                target_residue_list.push_back(buffer);
                getline(infile, buffer);
            }
            buffer.clear();
        }
    }

    Residue_Name_Conversion_Map conversion_map;
    if (isResidueConversionMap)
    {
        // Map to convert from NLN residue numbers in input file to HXB2 numbering via user provided map.
        std::string map_input_file = (workingDirectory + "/" + inputsDirectory + "/" + residue_conversion_map);
        conversion_map.ReadInputFile(map_input_file);
        //conversion_map.PrintMap();
    }

    // Load coordinates from all PDB files into one assembly, and create an assemblyvector too as stupid.
    MolecularModeling::Assembly assembly_of_pdbs(workingDirectory + "/" + inputsDirectory + "/" + template_pdb, gmml::InputFileType::PDB);
    assembly_of_pdbs.BuildStructureByDistance(4, 1.96);
    AssemblyVector vector_of_pdbs = loadPDBFilesFromDirectory(workingDirectory + "/" + snapshot_directory + "/", assembly_of_pdbs);
    AtomVector proteinBeads = beads::Add_Beads_To_Protein(assembly_of_pdbs);
    ResidueVector target_residues;
    ResidueVector all_residues = assembly_of_pdbs.GetAllResiduesOfAssembly();
    for(ResidueVector::iterator it1 = all_residues.begin(); it1 != all_residues.end(); ++it1)
    {
        MolecularModeling::Residue *current_residue = *it1;
        std::string id = current_residue->GetId();
        for(StringVector::iterator it2 = target_residue_list.begin(); it2 != target_residue_list.end(); ++it2)
        {
            std::string target_residue_number = *it2;
            std::string formatted_glycosite_number = "_" + target_residue_number + "_";
            if( id.compare(3, formatted_glycosite_number.size(), formatted_glycosite_number) == 0)
            {
                std::cout << "Bingo ya dingo: " << id << "\n";
                target_residues.push_back(current_residue);
            }
        }
    }

    MolecularModeling::Assembly receptors_ligand(workingDirectory + "/" + inputsDirectory + "/" + ligand_pdb, gmml::InputFileType::PDB);
    std::cout << "..." << std::endl;
    //CoordinateVector receptors_ligand_coordinates = receptors_ligand.GetAllCoordinates();
    CoordinateVector receptors_ligand_coordinates;
    for (auto &atom : receptors_ligand.GetAllAtomsOfAssembly())
    {
        if( (atom->GetName().compare("C1")==0) || (atom->GetName().compare("C3")==0) || (atom->GetName().compare("C5")==0) )
        {
            receptors_ligand_coordinates.push_back(atom->GetCoordinate());
        }
    }
    MolecularModeling::Assembly receptor(workingDirectory + "/" + inputsDirectory + "/" + receptor_pdb, gmml::InputFileType::PDB);
    std::cout << "...." << std::endl;
    AtomVector receptorBeads = beads::Add_Beads_To_Protein(receptor);
    std::cout << "....." << std::endl;
    CoordinateVector receptorCoordinates = receptor.GetAllCoordinates();
    std::cout << "......" << std::endl;


    Atom *receptorVector3DAtom1, *receptorVector3DAtom2;
    if(isReceptorVectorResidueList)
    {
        // Get atom pointers to CA atoms in the receptor (HA protein) that will be used for the HA oreintation vector
        receptorVector3DAtom1 = FindCAAtomOfResidueWithNumber(receptor.GetResidues(), receptor_vector_residue_list.at(0));
        receptorVector3DAtom2 = FindCAAtomOfResidueWithNumber(receptor.GetResidues(), receptor_vector_residue_list.at(1));
    }


    // For each target residue, go through each coordinate set (each snapshot) and get the coordinates of each atom.
    // Then, superimpose the ligand atoms to each snapshot with receptor and calculate overlap.
    bool reverse_strep_angle = false; // This is dumb, but reverse the angle every second residue? Real bad, not sure how else to do it.
    double overlap = 0.0;
    for(auto &target_residue : target_residues)
    {
        Atom *target_residue_C1_atom = target_residue->GetAtom("C1"); // Used in Streptavidin project to measure distance
        std::string converted_residue_number;
        AtomVector target_atoms = target_residue->GetAtoms();
        if (isResidueConversionMap)
        {
            converted_residue_number = conversion_map.ConvertViaMap(target_residue->GetNumber());
            std::cout << converted_residue_number << ": " << std::flush;
        }
        else
        {
            std::cout << target_residue->GetNumber() << ": " << std::flush;
            converted_residue_number = target_residue->GetNumber();
        }
        // When superimposing to glycans, the protein residue number is more useful to report:
//        MolecularModeling::Residue* aglycone = seek_alycone_residue_connected_to_residue(target_residue, &assembly_of_pdbs);
//        converted_residue_number = conversion_map.ConvertViaMap(aglycone->GetNumber());
//        std::cout << converted_residue_number << ": " << std::flush;

        unsigned int number_of_coordinate_sets = target_residue->GetAtoms().at(0)->GetCoordinates().size();
        std::cout << number_of_coordinate_sets << "\n";
        for (unsigned int i = 0; i < number_of_coordinate_sets; ++i)
        {
            //   MolecularModeling::Atom *current_atom = all_atoms[i];
            CoordinateVector target_atoms_coordinates;
            for(auto &atom : target_atoms)
            {
//            for(AtomVector::iterator it2 = target_atoms.begin(); it2!= target_atoms.end(); ++it2)
//            {
//                MolecularModeling::Atom *atom = *it2;
                if( (atom->GetName().compare("C1")==0) || (atom->GetName().compare("C3")==0) || (atom->GetName().compare("C5")==0) )
                {
                    target_atoms_coordinates.push_back(atom->GetCoordinates().at(i));
                }
            }
            // Now superimpose
            //std::cout << ".";
            //std::cout << receptors_ligand_coordinates.size() << "vs" << target_atoms_coordinates.size() << std::endl;
            gmml::Superimpose(receptors_ligand_coordinates, target_atoms_coordinates, receptorCoordinates);
            //std::cout << "..";

//            // WRITE OUT PDB FILES EVERY STEP
//            std::stringstream ss, ss1;
//           // ss << workingDirectory << "/receptor_" << i << "_" << aglycone->GetNumber() << ".pdb";
//            ss << workingDirectory << "/receptor_" << i << "_" << converted_residue_number << ".pdb";
//            PdbFileSpace::PdbFile *outputPdbFile = receptor.BuildPdbFileStructureFromAssembly(-1,0);
//            outputPdbFile->Write(ss.str());
//            ss1 << workingDirectory << "/ligand_" << i << "_" << converted_residue_number << ".pdb";
//            PdbFileSpace::PdbFile *outputPdbFile1 = receptors_ligand.BuildPdbFileStructureFromAssembly(-1, 0);
//            outputPdbFile1->Write(ss1.str());
            // END OF WRITE OUT PDB FILES EVERY STEP

            //std::cout << gmml::CalculateAtomicOverlaps(receptor.GetAllAtomsOfAssembly(), assembly_of_pdbs.GetAllAtomsOfAssembly()) << ", " << std::flush;
            overlap = beads::Calculate_bead_overlaps(receptorBeads, proteinBeads);
            std::cout << std::fixed << std::setprecision(3);
            std::cout << "Overlap: " << overlap << ", ";
            //std::cout << "." << std::flush;
            if(isReceptorVectorResidueList)
            {
                double angle;
                Eigen::Vector3d receptorVector3D(receptorVector3DAtom2->GetCoordinate()->GetX() - receptorVector3DAtom1->GetCoordinate()->GetX(),
                                                 receptorVector3DAtom2->GetCoordinate()->GetY() - receptorVector3DAtom1->GetCoordinate()->GetY(),
                                                 receptorVector3DAtom2->GetCoordinate()->GetZ() - receptorVector3DAtom1->GetCoordinate()->GetZ());
                if(reverse_strep_angle)
                    angle = GetAngleBetweenVectors(strepVector3DReversed, receptorVector3D, true);
                else
                    angle = GetAngleBetweenVectors(strepVector3D, receptorVector3D, true);

                // CALCULATE DISTANCE
                std::cout << "Angle: " << angle << ", ";
                std::cout << "GlycanToStrepMidlineDistance: " << GeometryTopology::calculateDistanceFromPointToLineBetweenTwoPoints(target_residue_C1_atom->GetCoordinates().at(i), streptavidinDisectingLinepoint1, streptavidinDisectingLinepoint2) << "\n";
           //     std::cout << "InterStrepDistance: " << receptorVector3DAtom2->GetDistanceToCoordinate(&streptavidinCenterPoint) << std::endl;


                std::cout << "InterStrepDistance: " << calculateInterStrepDistance(&streptavidinCenterPoint, receptorVector3DAtom1->GetCoordinate(), receptorVector3DAtom2->GetCoordinate());



                // Print out ATOM cards for testing in VMD
                std::cout << std::fixed << std::setprecision(3);
                std::cout << "ATOM      1 CA1  HA1     1    " << std::setw(8) << receptorVector3DAtom1->GetCoordinate()->GetX() << std::setw(8) << receptorVector3DAtom1->GetCoordinate()->GetY() << std::setw(8) << receptorVector3DAtom1->GetCoordinate()->GetZ() << "\n";
                std::cout << "ATOM      2 CA2  HA1     1    " << std::setw(8) << receptorVector3DAtom2->GetCoordinate()->GetX() << std::setw(8) << receptorVector3DAtom2->GetCoordinate()->GetY() << std::setw(8) << receptorVector3DAtom2->GetCoordinate()->GetZ() << "\n";
                std::cout << "ATOM      3 ST1  ST1     2    " << std::setw(8) << streptavidinDisectingLinepoint1.GetX() << std::setw(8) << streptavidinDisectingLinepoint1.GetY() << std::setw(8) << streptavidinDisectingLinepoint1.GetZ() << "\n";
                std::cout << "ATOM      4 ST2  ST1     2    " << std::setw(8) << streptavidinDisectingLinepoint2.GetX() << std::setw(8) << streptavidinDisectingLinepoint2.GetY() << std::setw(8) << streptavidinDisectingLinepoint2.GetZ() << "\n";
              //  std::cout << "ATOM      6 CA3  HA1     1    " << std::setw(8) << steve.GetX() << std::setw(8) << steve.GetY() << std::setw(8) << steve.GetZ() << "\n";

            }
        }
        reverse_strep_angle = !reverse_strep_angle; // Flips the bool from true to false for every target residue.
        std::cout << "\n";
    }

    //for(CoordinateVector::iterator it1 = )
    std::cout << "\n\nHEY OLY! Did you ensure that the Strep residues alternate from one side to the other? Cause you reversed the vector each time so that the angles would be the same. Don't forget that you derp.\n\n";
    return 0;
}

// Too tired for this crap. Just use the C1 it's the same.
//GeometryTopology::Coordinate getGeometricCenterOfCoordinates(GeometryTopology::CoordinateVector coordinates)
//{
//    double x_value, y_value, z_value = 0.0;
//    for (auto &coordinate : coordinates)
//    {
//        x_value += coordinate->GetX();
//        y_value += coordinate->GetY();
//        z_value += coordinate->GetZ();
//        GeometryTopology::Coordinate average_coordinate ( (x_value / coordinates.size()), (y_value / coordinates.size()), (z_value / coordinates.size()) );
//        return average_coordinate;
//    }
//}

//// Ignoring chain IDs. Rewrite or check other code for that.
//ResidueVector FindResiduesWithResidueNumber(ResidueVector residues, std::string query_residue_number)
//{
//    ResidueVector found_residues;
//    for(auto &residue : residues)
//    {
//        if (residue->GetNumber().compare(query_residue_number)==0)
//        {
//            found_residues.push_back(residue);
//        }
//    }
//    return found_residues;
//}

//Eigen::Vector3d calculateVector3DFromResidueSets(ResidueVector residueSet1, ResidueVector residueSet2)
//{
//    AtomVector atomSet1, atomSet2;
//    for(auto &residue : residueSet1)
//        atomSet1.push_back(residue->GetAtom("CA")); // Stop judging me I just want this to end.
//    for(auto &residue : residueSet2)
//        atomSet2.push_back(residue->GetAtom("CA"));

//    GeometryTopology::CoordinateVector
//}

double GetAngleBetweenVectors(Eigen::Vector3d v1, Eigen::Vector3d v2, bool returnDegrees)
{
    // angle = arccos ( v1.v2 / ||v1||||V2||) // "." is dot product. "||" is magnitude. norm() function gives magnitude
    if(returnDegrees)
        return ( (acos( v1.dot(v2) / (v1.norm() * v2.norm()))) * (180 / gmml::PI_RADIAN)); // calculate angle, then convert to degrees
    else
        return acos( v1.dot(v2) / (v1.norm() * v2.norm()));
}

Atom* FindCAAtomOfResidueWithNumber(ResidueVector residues, std::string query_residue_number)
{
    Atom* returning_atom_pointer = NULL;
    for(auto &residue : residues)
    {
        if (residue->GetNumber().compare(query_residue_number)==0)
        {
            returning_atom_pointer = residue->GetAtom("CA");
        }
    }
    return returning_atom_pointer;
}



AssemblyVector loadPDBFilesFromDirectory(std::string directory)
{
    //************************************************//
    // Load files from directory                      //
    //************************************************//

    std::cout << "loading files in directory: " << directory << std::endl;
    std::string filepath;
    AssemblyVector substrate_snapshots;

    DIR *dp; // A directory stream
    struct dirent *dirp; // Contains file serial number and name (char d_name[])
    struct stat filestat; // Contains info about file, such as device ID, user ID, access time etc

    dp = opendir( directory.c_str() ); //.c_str adds a null character to the end.
    if (dp == NULL)
    {
        std::cout << "Error(" << errno << ") opening " << directory << std::endl;
        return substrate_snapshots;
    }
    while ((dirp = readdir ( dp )))
    {
        filepath = directory + "/" + dirp->d_name;
        // If the file is a directory (or is in some way invalid) we'll skip it
        if (stat( filepath.c_str(), &filestat )) continue; // Is it a valid file?
        if (S_ISDIR( filestat.st_mode ))         continue; // Is it a directory?
        substrate_snapshots.push_back(new MolecularModeling::Assembly(filepath, gmml::InputFileType::PDB));
        std::cout << filepath << "\n";
    }
    closedir( dp );
    return substrate_snapshots;
}

AssemblyVector loadPDBFilesFromDirectory(std::string directory, MolecularModeling::Assembly combined_assembly)
{
    AssemblyVector pdb_files = loadPDBFilesFromDirectory(directory);
    std::cout << "Going through PDB files to copy coordinates into template.pdb\n";
    // Now go through and copy coordinates into one assembly?
    AtomVector master_atoms = combined_assembly.GetAllAtomsOfAssembly();
    for(AssemblyVector::iterator it1 = pdb_files.begin(); it1 != pdb_files.end(); ++it1)
    {
        MolecularModeling::Assembly* ass = *it1;
        AtomVector all_atoms = ass->GetAllAtomsOfAssembly();
        std::cout << "master.size() = " << master_atoms.size() << std::endl;
        std::cout << "all_atoms.size() = " << all_atoms.size() << std::endl;
        for (unsigned int i = 0; i < all_atoms.size(); ++i)
        {
            // push back coordinates of current atom onto corresponding atom's coordinate vector in combined_assembly
            MolecularModeling::Atom *current_atom = all_atoms[i];
            MolecularModeling::Atom *corresponding_master_atom = master_atoms[i];
            corresponding_master_atom->AddCoordinate(current_atom->GetCoordinate());
        }
    }
    return pdb_files;
}

MolecularModeling::Residue* seek_alycone_residue_connected_to_residue(MolecularModeling::Residue *residue, MolecularModeling::Assembly *assembly)
{
    bool found = false;
    bool *pfound = &found;
    MolecularModeling::Atom *start_atom = residue->GetAtoms().at(0); // Just get any atom to start from.
    MolecularModeling::Residue *aglycone = recursively_seek_alycone_residue(start_atom, pfound);
    //std::cout << "Finished with the aglycone: " << aglycone->GetId() << std::endl;

    //Clear the description; multiple calls to this function are possible.
    AtomVector temp = assembly->GetAllAtomsOfAssembly();
    for(AtomVector::iterator it = temp.begin(); it != temp.end(); ++it)
    {
        MolecularModeling::Atom *atom = *it;
        atom->SetDescription("");
    }
    return aglycone;
}

MolecularModeling::Residue* recursively_seek_alycone_residue(MolecularModeling::Atom *start_atom, bool *found)
{
    MolecularModeling::Residue *aglycone;
    start_atom->SetDescription("VisitedByAglyconSeeker"); // Hmm. Not sure this is a good idea, but it sure is handy.
   // std::cout << "Checking neighbors of " << start_atom->GetId() << std::endl;
    AtomVector neighbors = start_atom->GetNode()->GetNodeNeighbors();
    /*std::cout << "Neighbours are: " << std::endl;
    for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); ++it1)
    {
        std::cout << (*it1)->GetId() << std::endl;
    }
*/

    for(AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
    {
        // If haven't found aglycone yet and haven't visited this neighbor already
        if ( (*found == false) && ((*it)->GetDescription().compare("VisitedByAglyconSeeker")!=0) )
        {
            MolecularModeling::Atom *neighbor = *it;
            if(neighbor->GetResidue()->CheckIfProtein())
            {
                aglycone = neighbor->GetResidue();
                *found = true;
               // std::cout << "Found it!\n";
            }
            else
            {
                aglycone = recursively_seek_alycone_residue(neighbor, found);
            }
        }
    }
    return aglycone;
}

double calculateInterStrepDistance(GeometryTopology::Coordinate *streptavidinCenterPoint, GeometryTopology::Coordinate *influenzaHAPoint1, GeometryTopology::Coordinate *influenzaHAPoint2)
{
    /* HA1 and HA2 are coordinates in the center of the influenza HA tail and headgroup respectively.
 * ST3 is a coordinate at the center of the strep molecule.
 * Coordinates that start with a ? are ones that I want to calculate.
 * Once the HA is aligned to Strep via the above code, you can imagine an ST3mirror (ST3m) coordinate, that is a 120 degree rotation from ST3.
 * Why 120? The three HA binding sites are symetrically positioned around the HA headgroup, 120 degrees apart.
 * The purpose of the following code is to calculate the distance between the ST3m and ST3 (aka Interstrep Strep distance on surface required for binding like this).
 * HA3 is a coordinate that is in line with the HA1->HA2 vector and in the same plane as ST3 and ST3m, such that ST3-HA3-ST3m forms a 120 degree angle and an isoceles triangle.
 * I first calculate the position of HA3:
 * I get the HA2-ST3 distance. The HA1-HA2-ST3 angle, and calculate HA2-HA3 from that, as HA2-HA3-ST3 angle is 90. Math yo.
 * I then calculate HA3 position using the normalized HA1-HA2 vector and the HA2-HA3 distance.
 * To calculate ST3m - ST3 distance, I use the HA3 coordinate and the fact that it is an isoceles triangle with a ST3-HA3-ST3m angle of 120. I can calculate HA3-ST3 distance too.
 *
 *                 ?ST3mirror
 *
 *  HA1--------HA2 ?HA3
 *
 *                  ST3 // Streptavidin center point
 *
 */
    Eigen::Vector3d influenzaHAVector3D(influenzaHAPoint2->GetX() - influenzaHAPoint1->GetX(),
                                        influenzaHAPoint2->GetY() - influenzaHAPoint1->GetY(),
                                        influenzaHAPoint2->GetZ() - influenzaHAPoint1->GetZ());

    Eigen::Vector3d ha2ToStrepCenter(streptavidinCenterPoint->GetX() - influenzaHAPoint2->GetX(),
                                     streptavidinCenterPoint->GetY() - influenzaHAPoint2->GetY(),
                                     streptavidinCenterPoint->GetZ() - influenzaHAPoint2->GetZ());

    double angleHA2StrepCenter = GetAngleBetweenVectors(influenzaHAVector3D, ha2ToStrepCenter);
    std::cout << "Outer HA1 HA2 StrepCenter Angle: " << GetAngleBetweenVectors(influenzaHAVector3D, ha2ToStrepCenter, true) << "\n";

    double hyp1 = influenzaHAPoint2->Distance(streptavidinCenterPoint);
    std::cout << "Cos result of " << GetAngleBetweenVectors(influenzaHAVector3D, ha2ToStrepCenter, true) << " is " << cos(angleHA2StrepCenter) << "\n";
    double adj1 = hyp1 * cos(angleHA2StrepCenter);
    std::cout << "adj1: " << adj1 << " based on hyp of " << hyp1 << "\n";

    double lengthRatio = ( (adj1 + influenzaHAVector3D.norm()) / influenzaHAVector3D.norm());
    GeometryTopology::Coordinate ha3(((1-lengthRatio) * influenzaHAPoint1->GetX()) + lengthRatio * influenzaHAPoint2->GetX(),
                                        ((1-lengthRatio) * influenzaHAPoint1->GetY()) + lengthRatio * influenzaHAPoint2->GetY(),
                                        ((1-lengthRatio) * influenzaHAPoint1->GetZ()) + lengthRatio * influenzaHAPoint2->GetZ());

    double hyp2 = ha3.Distance(streptavidinCenterPoint);
    double opp2 = hyp2 * sin(1.0472); // Cos of 60 degrees.
    std::cout << "Distance between Strep Centers if finally: " << (opp2 * 2) << ". Based on hyp2: " << hyp2 << "\n";
    std::cout << "ATOM      6 CA4  HA1     1    " << std::setw(8) << ha3.GetX() << std::setw(8) << ha3.GetY() << std::setw(8) << ha3.GetZ() << "\n";
    return (opp2 * 2);

}


