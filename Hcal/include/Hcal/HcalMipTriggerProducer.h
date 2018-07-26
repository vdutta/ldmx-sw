/**
 * @file HcalMipTriggerProducer.h
 * @brief
 * @author
 */

#ifndef HCAL_HCALMIPTRIGGERPRODUCER_H
#define HCAL_HCALMIPTRIGGERPRODUCER_H

//LDMX Framework
#include "Event/Event.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file

namespace ldmx {
    
    /**
     * @class HcalMipTriggerProducer
     * @brief 
     */
    class HcalMipTriggerProducer : public ldmx::Producer {
        public:

            HcalMipTriggerProducer(const std::string& name, ldmx::Process& process) : ldmx::Producer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void produce(const ldmx::Event& event);

            virtual void onFileOpen();

            virtual void onFileClose();

            virtual void onProcessStart(); 

            virtual void onProcessEnd();

        private:

    };
}

#endif /* HCAL_HCALMIPTRIGGERPRODUCER_H */
