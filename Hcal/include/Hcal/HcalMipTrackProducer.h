/**
 * @file HcalMipTrackProducer.h
 * @brief Header file for HcalMipTrackProducer class
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_HCALMIPTRACKPRODUCER_H
#define HCAL_HCALMIPTRACKPRODUCER_H

//LDMX Framework
#include "Event/Event.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file

namespace ldmx {
    
    /**
     * @class HcalMipTrackProducer
     * @brief Producer that reconstructs MIP tracks through the hcal.
     */
    class HcalMipTrackProducer : public ldmx::Producer {
        public:

            HcalMipTrackProducer(const std::string& name, ldmx::Process& process) : ldmx::Producer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void produce(ldmx::Event& event);

            virtual void onFileOpen() { }

            virtual void onFileClose() { }

            virtual void onProcessStart() { } 

            virtual void onProcessEnd() { }

        private:

    };
}

#endif /* HCAL_HCALMIPTRACKPRODUCER_H */
