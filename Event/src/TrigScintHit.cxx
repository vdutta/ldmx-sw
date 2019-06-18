#include "Event/TrigScintHit.h"

// STL
#include <iostream>

ClassImp(ldmx::TrigScintHit)

namespace ldmx {

    void TrigScintHit::Clear(Option_t *option) {
        HcalHit::Clear();
        setPE(0);
    }

    void TrigScintHit::Print(Option_t *option) const {
        std::cout << "TrigScintHit { " << "id: " << std::hex << getID() << std::dec
                << ",  energy: " << getEnergy() << "MeV, time: " << getTime()
                << "ns, amplitude: " << getAmplitude() << ", pe: " << getPE() << "}" << std::endl;
    }

}
