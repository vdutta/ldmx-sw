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
 * @app ldmx-charm-gdml
 * "Charms" gdml files into being a single TGeo-friendly gdml file.
 *
 * Run this exactly the same way you would run ldmx-sim but without a
 * steering macro. It pulls the geometry description from the gdml
 * files the same way ldmx-sim does, so you need to make sure that
 * all of the files you want are in your current directory.
 * This utility assumes the parent gdml file is titled 'detector.gdml'.
 */
int main(int nargs, const char** )  {
    
    //Validate inputs
    if ( nargs > 1 ) {
        //Some input was given, assume user is requesting help
        printf("Usage: ldmx-charm-gdml\n" );
        printf("    No inputs are needed, but the necessary gdml files need to be in the\n" );
        printf("    current working directory. It is suggested that you sym-link the    \n" );
        printf("    detector gdml files to the directory you are in using 'ln -s'.      \n" );
        return 1;
    }

    //Import geometry to Geant4

    //parse detector into Geant4
    G4GDMLParser read_parser;
    read_parser.Read("detector.gdml");

    //get world volume from parser
    G4VPhysicalVolume* worldVol = read_parser.GetWorldVolume();

    //Export geometry to new file
    //  Separate parser to force all geometries written to one file
    G4GDMLParser write_parser;
    write_parser.Write( 
            "TGeofriend_detector.gdml" , 
            worldVol 
            );

    return 0;
}
