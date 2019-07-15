#include "SimApplication/G4TriggerPadHit.h"

// Geant4
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Point3D.hh"
#include "G4Circle.hh"

namespace ldmx {

    G4Allocator<G4TriggerPadHit> G4TriggerPadHitAllocator;

    void G4TriggerPadHit::Draw() {

        G4VVisManager* visManager = G4VVisManager::GetConcreteInstance();

        if (visManager) {

            G4Point3D p3D = G4Point3D(this->position_);
            G4Circle chit(p3D);
            chit.SetScreenDiameter(3.0);
            chit.SetFillStyle(G4Circle::filled);

            G4Colour col(1.0, 0.0, 1.0);

            chit.SetVisAttributes(G4VisAttributes(col));
            visManager->Draw(chit);
        }
    }

    void G4TriggerPadHit::Print() {
        print(std::cout);
    }

    std::ostream& G4TriggerPadHit::print(std::ostream& os) {
        os << "G4TriggerPadHit { " "edep: " << this->edep_
                << ", " "position: (" << this->position_[0] << ", " << this->position_[1] << ", " << this->position_[2] << "), "
                <<  "padID: " << this->padID_ << ", "
                <<  "layerID: " << this->layerID_ << ", "
                <<  "strip: " << this->stripID_ << ", "
                << "momentum: (" << this->momentum_[0] << ", " << this->momentum_[1] << ", " << this->momentum_[2] << "), "
	  //     << "pathLength: " << this->pathLength_ << ", " 
	        << "time: " << this->time_ << " [ns] }" << std::endl;
        return os;
    }

}

