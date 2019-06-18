#include "Event/SimTriggerPadHit.h"

ClassImp(ldmx::SimTriggerPadHit)

namespace ldmx {

    SimTriggerPadHit::SimTriggerPadHit() : TObject() {
    }

    SimTriggerPadHit::~SimTriggerPadHit() {
        Clear();
    }

    void SimTriggerPadHit::Print(Option_t *option) const {
        std::cout << "SimTriggerPadHit { " << "id: " << id_ << ", " <<
                "stripID: " << stripID_ << ", " <<
                "layerID: " << layerID_ << ", " <<
                "padID: " << padID_ << ", " <<
                "position: ( " << x_ << ", " << y_ << ", " << z_ << " ), " <<
                "edep: " << edep_ << ", " <<
                "time: " << time_ << ", " <<
                "momentum: ( " << px_ << ", " << py_ << ", " << pz_ << " )" <<
                " }" << std::endl;
    }

    void SimTriggerPadHit::Clear(Option_t*) {
        TObject::Clear();

        id_ = 0;
        stripID_ = 0;
        layerID_ = 0;
        padID_ = 0;
        edep_ = 0;
        time_ = 0;
        px_ = 0;
        py_ = 0;
        pz_ = 0;
        x_ = 0;
        y_ = 0;
        z_ = 0;
        energy_ = 0;
        pathLength_ = 0;
        trackID_ = -1;
        pdgID_ = 0;

        simParticle_ = nullptr;
    }

    SimParticle* SimTriggerPadHit::getSimParticle() const {
        return static_cast<SimParticle*>(simParticle_.GetObject());
    }

    void SimTriggerPadHit::setPosition(const float x, const float y, const float z) {
        this->x_ = x;
        this->y_ = y;
        this->z_ = z;
    }

    void SimTriggerPadHit::setMomentum(const float px, const float py, const float pz) {
        this->px_ = px;
        this->py_ = py;
        this->pz_ = pz;
    }
}
