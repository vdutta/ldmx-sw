/**
 * @file TriggerPadSD.h
 * @brief Class defining a trigger scintillator pad sensitive detector
 * @author Jeremy McCormick, SLAC National Accelerator Laboratory; adapted from Hcal by Lene Kristian Bryngemark, Stanford 
 */

#ifndef SIMAPPLICATION_TRIGGERPADSD_H_
#define SIMAPPLICATION_TRIGGERPADSD_H_

// LDMX
#include "SimApplication/CalorimeterSD.h"
//#include "DetDescr/DetectorID.h"
#include "DetDescr/TriggerPadID.h"
//#include "SimApplication/G4TrackerHit.h"
#include "SimApplication/G4TriggerPadHit.h"

namespace ldmx {

    /**
     * @class TriggerPadSD
     * @brief TriggerPad sensitive detector
     *
     * @note
     * This class sets hit coordinates, pdgid association, and detector id.
     *
     * @todo Add actual custom hit processing for Trigger scintillator pad detector.
     */
  class TriggerPadSD : public CalorimeterSD { // public G4VSensitiveDetector { 

        public:

            TriggerPadSD(G4String name, G4String theCollectionName, int subdet, DetectorID* detID = new TriggerPadID);

            virtual ~TriggerPadSD();

            G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
            
        private:
          double birksc1_;
          double birksc2_;    

	  /**                                                                                                                                                                               
	   * The output hits collection of G4TrackerHits.                                                                                                                                   
	   */
	  G4TriggerPadHitsCollection* hitsCollection_;

	  /**                                                                                                                                                                               
	   * The subdetector ID.                                                                                                                                                            
	   */
	  int subdetID_;

	  /**                                                                                                                                                                               
	   * The detector ID.                                                                                                                                                               
	   */
	  DetectorID* detID_{new TriggerPadID};


	  /**
	   * Initialize the sensitive detector.
	   * @param hcEvent The hits collections of the event.
	   */
	  void Initialize(G4HCofThisEvent* hcEvent);

	  /**                                                                                                                                                                               
	   * End of event hook.                                                                                                                                                             
	   * @param hcEvent The hits collections of the event.                                                                                                                              
	   */
	  void EndOfEvent(G4HCofThisEvent* hcEvent);


    };

}

#endif
