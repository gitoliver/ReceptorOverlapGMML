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
AssemblyVector loadPDBFilesFromDirectory(std::string directory);
AssemblyVector loadPDBFilesFromDirectory(std::string directory, MolecularModeling::Assembly combined_assembly);
MolecularModeling::Residue* seek_alycone_residue_connected_to_residue(MolecularModeling::Residue *residue, MolecularModeling::Assembly *assembly);
MolecularModeling::Residue* recursively_seek_alycone_residue(MolecularModeling::Atom *start_atom, bool *found);
//GeometryTopology::Coordinate getGeometricCenterOfCoordinates(CoordinateVector coordinates);

int main(int argc, char *argv[])
{



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

    GeometryTopology::Coordinate receptorDisectingLinepoint1(40.1895, 21.084, 7.041);
    GeometryTopology::Coordinate receptorDisectingLinepoint2(6.096, 21.8635, -6.911);
    GeometryTopology::Coordinate jola(1.0, -15.0, 1.0);

    std::cout << "Distance from jola to line is " << GeometryTopology::calculateDistanceFromPointToLineBetweenTwoPoints(jola, receptorDisectingLinepoint1, receptorDisectingLinepoint2) << "\n";

    typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;

    std::string workingDirectory = Find_Program_Working_Directory();
    std::string inputsDirectory = argv[1];
    std::string inputFile = argv[2];
    //************************************************//
    // Reading input file                             //
    //************************************************//

    StringVector target_residue_list;
    std::string buffer, snapshot_directory, template_pdb, ligand_pdb, receptor_pdb, residue_conversion_map;
    bool isResidueConversionMap = false;
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


    // For each target residue, go through each coordinate set (each snapshot) and get the coordinates of each atom.
    // Then, superimpose the ligand atoms to each snapshot with receptor and calculate overlap.
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
            std::cout << GeometryTopology::calculateDistanceFromPointToLineBetweenTwoPoints(target_residue_C1_atom->GetCoordinates().at(i), receptorDisectingLinepoint1, receptorDisectingLinepoint2) << ", " << std::flush;
            //std::cout << "..";
//            std::stringstream ss, ss1;
//           // ss << workingDirectory << "/receptor_" << i << "_" << aglycone->GetNumber() << ".pdb";
//            ss << workingDirectory << "/receptor_" << i << "_" << converted_residue_number << ".pdb";
//            PdbFileSpace::PdbFile *outputPdbFile = receptor.BuildPdbFileStructureFromAssembly(-1,0);
//            outputPdbFile->Write(ss.str());
//            ss1 << workingDirectory << "/ligand_" << i << "_" << converted_residue_number << ".pdb";
//            PdbFileSpace::PdbFile *outputPdbFile1 = receptors_ligand.BuildPdbFileStructureFromAssembly(-1, 0);
//            outputPdbFile1->Write(ss1.str());
            //std::cout << gmml::CalculateAtomicOverlaps(receptor.GetAllAtomsOfAssembly(), assembly_of_pdbs.GetAllAtomsOfAssembly()) << ", " << std::flush;
            overlap = beads::Calculate_bead_overlaps(receptorBeads, proteinBeads);
            std::cout << overlap << ", " << "\n";
//            if (overlap <= 0.5)
//            {
//                // Get shapes with good HA orientation
//                // Calculate strep intersection line
//                //
//                // Calculate minimum distance from point to line
//                //GeometryTopology::calculateDistanceFromPointToLineBetweenTwoPoints(queryPoint, linePointA, linePointB);
//            }
        }
        std::cout << "\n";
    }

    //for(CoordinateVector::iterator it1 = )

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




