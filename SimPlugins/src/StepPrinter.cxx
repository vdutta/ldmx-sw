/**
 * @file StepPrinter.cxx
 * @class StepPrinter
 * @brief Sim plugin that prints the details of a step taken by a particle.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#include "SimPlugins/StepPrinter.h"

//------------//
//   Geant4   //
//------------//
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

//----------------//
//   C++ StdLib   //
//----------------//
#include <iostream>

SIM_PLUGIN(ldmx, StepPrinter)

namespace ldmx { 

    G4ClassificationOfNewTrack StepPrinter::stackingClassifyNewTrack(
            const G4Track* track, 
            const G4ClassificationOfNewTrack& currentTrackClass) {

        // Get the particle type.
        G4String particleName = track->GetParticleDefinition()->GetParticleName();
       
        // get the PDGID of the track.
        G4int pdgID = track->GetParticleDefinition()->GetPDGEncoding();

        std::string process{""}; 
        if (track->GetParentID() != 0) process = track->GetCreatorProcess()->GetProcessName();

        std::cout << "*******************************\n" 
                  << "*   Pushing to stack\n"  
                  << "********************************\n"
                  << "\tParticle " << particleName << " ( PDG ID: " << pdgID << " ) : " << "\n"
                  << "\tTrack ID: " << track->GetTrackID() << "\n"
                  << "\tCreated via " << process << "\n" 
                  << "\tKinetic Energy: " << track->GetKineticEnergy() << "\n"
                  << "********************************\n"
                  << "********************************"
                  << std::endl;


        return currentTrackClass; 
    }

    void StepPrinter::stepping(const G4Step* step) {
       
        // Get the track associated with this step.
        G4Track* track = step->GetTrack();

        // get the PDGID of the track.
        G4int pdgID = track->GetParticleDefinition()->GetPDGEncoding();

        // Get the volume the particle is in.
        G4String volumeName = track->GetVolume()->GetName();

        // Get the particle type.
        G4String particleName = track->GetParticleDefinition()->GetParticleName();

        // Get the kinetic energy of the particle.
        double ke = track->GetKineticEnergy();

        std::string process{""}; 
        if (track->GetParentID() != 0) process = track->GetCreatorProcess()->GetProcessName();

        std::cout << "*******************************\n" 
                  << "*   Step " << track->GetCurrentStepNumber() << "\n"
                  << "********************************\n"
                  << "[ StepPrinter ]: "    << "\n" 
                  << "\tTotal energy of "   << particleName 
                  << " ( PDG ID: " << pdgID << " ) : " << ke << "\n"
                  << "\tParticle currently in " << volumeName << "\n"
                  << "\tTrack ID: " << track->GetTrackID() << "\n" 
                  << "\tStep #: " << track->GetCurrentStepNumber() << "\n"
                  << "\tSecondaries: " << step->GetSecondaryInCurrentStep()->size() << "\n"
                  << "\tCreated via " << process << "\n" 
                  << "\tPost step process: " << step->GetPostStepPoint()->GetStepStatus() << "\n" 
                  << "********************************\n"
                  << "********************************"
                  << std::endl;
    }
}
