/*******************************************************************
 * @file ldmx_charm_gdml.cxx
 * @author Tom Eichlersmith, University of Minnesota
 *
 * This file contains the utility that makes ROOT-friendly gdml
 * ("charms" the gdml) from more complicated gdml files.
 * This utility uses Geant4 to parse the input detector
 * geometry and the re-export it to a file containing all of the
 * detector geometry without the more complicated gdml bits (e.g. loops
 * and replica volumes).
 */

///////////////
// Geant4
//////////////
#include "G4GDMLParser.hh"
#include "G4VPhysicalVolume.hh"

////////////////
// STL
///////////////
#include <iostream> //for printf

/********************************************************************
 * @function print_help
 * Prints help message for ldmx-charm-gdml
 */
void print_help_ldmx_charm_gdml();

/********************************************************************
 * @app ldmx-charm-gdml
 * "Charms" gdml files into being a single TGeo-friendly gdml file.
 *
 * Run this exactly the same way you would run ldmx-sim but without a
 * steering macro. It pulls the geometry description from the gdml
 * files the same way ldmx-sim does, so you need to make sure that
 * all of the files you want are in your current directory.
 * This utility assumes the parent gdml file is titled 'detector.gdml'.
 *
 * You can optionally input a gdml file to charm but make sure its
 * dependencies are in your current working directory (if there are
 * any).
 */
int main(int nargs, const char** argv )  {
    
    //Validate inputs
    std::string gdml_file_name = "detector.gdml";
    if ( nargs > 2 ) {
        //Extra input was given, assume user is requesting help
        print_help_ldmx_charm_gdml();
        return 1;
    } else if ( nargs == 2 ) {
        //gdml-file-name argument is given
        gdml_file_name = argv[1];

        std::size_t help_pos = gdml_file_name.find("help");
        //check if asking for help
        if ( 
                gdml_file_name == "h" or
                gdml_file_name == "-h" or
                help_pos != std::string::npos
           ) {
            print_help_ldmx_charm_gdml();
            return 0;
        }//print help and exit
    }//input was given

    /////////////////
    //Import geometry to Geant4

    //parse detector into Geant4
    G4GDMLParser read_parser;
    read_parser.Read( gdml_file_name.c_str() );

    //get world volume from parser
    G4VPhysicalVolume* worldVol = read_parser.GetWorldVolume();

    //////////////////////////
    //Export geometry to new file

    //output file same name but with prefix
    std::string output_prefix = "TGeoFriend_";
    std::string output_name = output_prefix + gdml_file_name;
    //separate parser to force all geometries written to one file
    G4GDMLParser write_parser;
    write_parser.Write( output_name.c_str() , worldVol );

    return 0;
}

void print_help_ldmx_charm_gdml() {
    printf("Usage: ldmx-charm-gdml [gdml-file-path]\n" );
    printf(" [gdml-file-name]: name of gdml file you wish to charm (OPTIONAL - default is 'detector.gdml')\n" );
    printf("    MUST be in current working directory and have '.gdml' extension\n" );
    printf("\n");
    printf("    No inputs are required, but the necessary gdml files need to be in the\n" );
    printf("    current working directory. It is suggested that you sym-link the    \n" );
    printf("    detector gdml files to the directory you are in using 'ln -s'.      \n" );
    return;
}
