/**
 * @file LHEPrimaryGenerator.h
 * @brief Class for generating a Geant4 event from LHE event data
 * @author Jeremy McCormick, SLAC National Accelerator Laboratory
 */

#ifndef SIMAPPLICATION_LHEPRIMARYGENERATOR_H_
#define SIMAPPLICATION_LHEPRIMARYGENERATOR_H_

// Geant4
#include "G4RunManager.hh"
#include "G4VPrimaryGenerator.hh"

// LDMX
//#include "SimApplication/LHEReader.h"
#include "SimApplication/LHEF.h"

// STL
#include <string> //lhe file name
#include <sstream> //parsing #vertex flag

namespace ldmx {

    /**
     * @class LHEPrimaryGenerator
     * @brief Generates a Geant4 event from an LHEEvent
     */
    class LHEPrimaryGenerator : public G4VPrimaryGenerator {

        public:

            /**
             * Class constructor.
             * @param reader The LHE reader with the event data.
             */
            LHEPrimaryGenerator( std::string lhe_file_name );

            /**
             * Class destructor.
             */
            virtual ~LHEPrimaryGenerator();

            /**
             * Generate vertices in the Geant4 event.
             * @param anEvent The Geant4 event.
             */
            void GeneratePrimaryVertex(G4Event* anEvent);

        private:

            /**
             * The LHE reader with the event data.
             */
            LHEF::Reader reader_;
    };

}

#endif
