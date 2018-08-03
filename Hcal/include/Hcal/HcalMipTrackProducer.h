/**
 * @file HcalMipTrackProducer.h
 * @brief Header file for HcalMipTrackProducer class
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_HCALMIPTRACKPRODUCER_H
#define HCAL_HCALMIPTRACKPRODUCER_H

//STL
#include <string> //names of collections

//LDMX Framework
#include "Event/Event.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "Event/EventConstants.h" //for HcalMipTrack string
#include "Event/HcalMipTrack.h" //mip track container

namespace ldmx {
    
    /**
     * @class HcalMipTrackProducer
     * @brief Producer that reconstructs MIP tracks through the hcal.
     */
    class HcalMipTrackProducer : public ldmx::Producer {
        public:

            HcalMipTrackProducer(const std::string& name, ldmx::Process& process) : ldmx::Producer(name, process) { }

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void produce(ldmx::Event& event);

            virtual void onFileOpen() { }

            virtual void onFileClose() { }

            virtual void onProcessStart() { } 

            virtual void onProcessEnd() { }

        private:
            
            /** Name of collection of HcalHits */
            std::string hcalHitCollName_;

            /** Name of pass to get HcalHits collection from */
            std::string hcalHitPassName_;

            /** Container for reconstructed tracks */
            TClonesArray *hcalMipTracks_;

            /** Collection name for event bus */
            std::string hcalMipTracksCollName_;

            /** Radius of cylinder around track [mm] */
            double trackRadius_;

            /** Minimum number of PE to be considered not noise */
            double minPE_;

            /** Maximum energy of a cluster to be considered a mip */
            double maxEnergy_;

            /** Minimum number of clusters in track to be considered for ranking */
            int minNumClusters_;

            
    };
}

#endif /* HCAL_HCALMIPTRACKPRODUCER_H */
