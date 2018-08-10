/**
 * @file MipChecker.h
 * @brief Header file for class MipChecker
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_MIPCHECKER_H
#define HCAL_MIPCHECKER_H

#include "TClonesArray.h"

//LDMX Framework
#include "Event/Event.h"
#include "Event/TriggerResult.h"
#include "Event/HcalMipTrack.h"
#include "Event/SimCalorimeterHit.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file

namespace ldmx {
    
    /**
     * @class MipChecker
     * @brief Class to compare mip track and trigger producers to actual sim particles.
     */
    class MipChecker : public ldmx::Analyzer {
        public:

            MipChecker(const std::string& name, ldmx::Process& process) : ldmx::Analyzer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void analyze(const ldmx::Event& event);

            virtual void onFileOpen() { }

            virtual void onFileClose() { }

            virtual void onProcessStart() { }

            virtual void onProcessEnd();

        private:
            
            /** Counters for trigger */
            unsigned int numFalsePos_;
            unsigned int numTruePos_;
            unsigned int numFalseNeg_;
            unsigned int numTrueNeg_;

            /** Counters for track */
            unsigned int numWrong_;
            unsigned int numRight_;

            /** Number simulated events */
            unsigned int numEvents_;
    };
}

#endif /* HCAL_MIPCHECKER_H */
