#include "SimApplication/TriggerPadSD.h"

// STL
#include <iostream>

// Geant4
#include "G4Box.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SDManager.hh"
#include "G4Geantino.hh"
#include "G4ChargedGeantino.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

namespace ldmx {

    TriggerPadSD::TriggerPadSD(G4String name, G4String theCollectionName, int subdetID, DetectorID* detID) :
            CalorimeterSD(name, theCollectionName, subdetID, detID),
            birksc1_(1.29e-2),
            birksc2_(9.59e-6)
    {
    }

    TriggerPadSD::~TriggerPadSD() {}

    G4bool TriggerPadSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist) {

        // Determine if current particle of this step is a Geantino.
        G4ParticleDefinition* pdef = aStep->GetTrack()->GetDefinition();
        bool isGeantino = false;
        if (pdef == G4Geantino::Definition() || pdef == G4ChargedGeantino::Definition()) {
            isGeantino = true;
        }

        // Get the edep from the step.
        G4double edep = aStep->GetTotalEnergyDeposit();

        // Skip steps with no energy dep which come from non-Geantino particles.
        if (edep == 0.0 && !isGeantino) {
            if (verboseLevel > 2) {
                std::cout << "TriggerPadSD skipping step with zero edep." << std::endl << std::endl;
            }
            return false;
        }

        //---------------------------------------------------------------------------------------------------
        //                Birks' Law 
        //                ===========
        //    
        //                      In the case of Scintillator as active medium, we can
        //              describe the quenching effects with the Birks' law,
        //              using the expression and the coefficients taken from
        //              the paper NIM 80 (1970) 239-244 for the organic
        //              scintillator NE-102:
        //                                     S*dE/dr
        //                    dL/dr = -----------------------------------
        //                               1 + C1*(dE/dr)
        //              with:
        //                    S=1
        //                    C1 = 1.29 x 10^-2  g*cm^-2*MeV^-1
        //                    C2 = 9.59 x 10^-6  g^2*cm^-4*MeV^-2
        //              These are the same values used by ATLAS TileCal
        //              and CMS HCAL (and also the default in Geant3).
        //              You can try different values for these parameters,
        //              to have an idea on the uncertainties due to them,
        //              by uncommenting one of the lines below.
        //              To get the "dE/dr" that appears in the formula,
        //              which has the dimensions
        //                    [ dE/dr ] = MeV * cm^2 / g
        //              we have to divide the energy deposit in MeV by the
        //              product of the step length (in cm) and the density
        //              of the scintillator:

        G4double birksFactor(1.0);
        G4double stepLength = aStep->GetStepLength() / CLHEP::cm;    
        //Do not apply Birks for gamma deposits!
        if (stepLength > 1.0e-6) //Check, cut if necessary.
        {
            G4double rho = aStep->GetPreStepPoint()->GetMaterial()->GetDensity() / (CLHEP::g / CLHEP::cm3);
            G4double dedx = edep / (rho * stepLength); //[MeV*cm^2/g]
            birksFactor = 1.0 / (1.0 + birksc1_ * dedx + birksc2_ * dedx * dedx);
            if (aStep->GetTrack()->GetDefinition() == G4Gamma::GammaDefinition()) birksFactor = 1.0;
            if (aStep->GetTrack()->GetDefinition() == G4Neutron::NeutronDefinition()) birksFactor = 1.0;  
        }


        // Create a new hit. 
        G4TriggerPadHit* hit = new G4TriggerPadHit();

        // Set the edep.
        hit->setEdep(edep*birksFactor);
 
        // Get the scintillator solid box
	//        G4Box* scint = static_cast<G4Box*>(aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetSolid());

        // Set the step mid-point as the hit position.
        G4StepPoint* prePoint = aStep->GetPreStepPoint();
        G4StepPoint* postPoint = aStep->GetPostStepPoint();
        G4ThreeVector position = 0.5 * (prePoint->GetPosition() + postPoint->GetPosition());
	//        G4ThreeVector localPosition = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(position);        
        hit->setPosition(position[0], position[1], position[2]);
        
        // Set the global time.
        hit->setTime(aStep->GetTrack()->GetGlobalTime());

	G4ThreeVector p = postPoint->GetMomentum();
	hit->setMomentum(p.x(), p.y(), p.z());


        // Create the ID for the hit.
        int copyNum = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetVolume(layerDepth_)->GetCopyNo();
        int padID = copyNum / 100;
	//alternate numbering means layers are all odd/even copy numbers
        int layer = copyNum % 2;
        int stripID = copyNum % 100;

        detID_->setFieldValue(1, padID);
        detID_->setFieldValue(2, layer);
        detID_->setFieldValue(3, stripID);
        hit->setID(detID_->pack());

        // Set the track ID on the hit.
        hit->setTrackID(aStep->GetTrack()->GetTrackID());

        // Set the PDG code from the track.
        hit->setPdgID(aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding());


        if (this->verboseLevel > 2) {
            std::cout << "Created new SimHit in detector " << this->GetName()
                      << " subdet ID <" << subdet_ << ">, layer <" << layer << "> and padID <" << padID << ">, copynum <" << copyNum << ">"
                      << std::endl;
            hit->Print();
            std::cout << std::endl;
        }

        // Insert the hit into the hits collection.
        hitsCollection_->insert(hit);

        return true;
    }
  void TriggerPadSD::Initialize(G4HCofThisEvent* hce) {

    // Setup hits collection and the HC ID.                                                                                                                                               
    hitsCollection_ = new G4TriggerPadHitsCollection(SensitiveDetectorName, collectionName[0]);
    int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(hcID, hitsCollection_);
  }
  
}

