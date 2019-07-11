#include "../includes/residue_name_conversion_map.h"
Residue_Name_Conversion_Map::Residue_Name_Conversion_Map()
{
}

Residue_Name_Conversion_Map::Residue_Name_Conversion_Map(std::string input_file)
{
    ReadInputFile(input_file);
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

std::string Residue_Name_Conversion_Map::ConvertViaMap(std::string residue_number)
{
    for(PairVector::iterator it = map_.begin(); it != map_.end(); ++it)
    {
       //std::cout << "Comparing " << residue_number << " with " << (*it).first << std::endl;
       if (residue_number.compare((*it).first)==0)
       {
            return (*it).second;
       }
    }
    std::cout << "Error, did not find residue number in map\n";
    return residue_number;
}


//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

void Residue_Name_Conversion_Map::ReadInputFile(std::string input_file)
{
    std::ifstream infile (input_file);
    if (!infile)
    {
        std::cerr << "Uh oh, input file could not be opened for reading!" << std::endl;
        std::exit(1);
    }
    while (infile) // While there's still stuff left to read
    {
        std::string strInput;
        getline(infile, strInput);
        //std::cout << "Line is " << strInput << std::endl;
        StringVector splitLine = split(strInput, ',');
        std::string second;
        for(StringVector::iterator it = splitLine.begin(); it != splitLine.end(); ++it)
        {
            if (it == splitLine.begin())
            {
                second = (*it);
            }
            else
            {
                map_.push_back(make_pair((*it), second));
            }
        }
    }
    return;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

void Residue_Name_Conversion_Map::PrintMap()
{
    for(PairVector::iterator it = map_.begin(); it != map_.end(); ++it)
    {
        std::cout << it->first << "=" << it->second << " " << std::endl;
    }
    return;
}

