/**
 * @file HcalTrackAnalyzer.h
 * @brief Class implementation of analyzer to study track differences in HCal
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef HCAL_HCALTRACKANALYZER_H
#define HCAL_HCALTRACKANALYZER_H

//Standard Libraries
#include <iostream> //Checks to std::cout for development purposes
#include <string> //names of histograms

//ROOT
#include "TH1.h" //One Dimensional Histograms
#include "TClonesArray.h" //Array to get list of hits from Event

//LDMX Framework
#include "Framework/EventProcessor.h" //Needed to declare analyzer
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "Event/Event.h" //Study event by event
#include "Event/HcalMipTrack.h" //Study hcal tracks
#include "Event/SimParticle.h" //check actual number of muons
#include "Event/HcalHit.h" //check if enough hcal hits to reconstruct

namespace ldmx {
    
    /**
     * @class HcalTrackAnalyzer
     * @brief ldmx::Analyzer that constructs histograms studying how tracks in the Hcal behave differently.
     */
    class HcalTrackAnalyzer : public ldmx::Analyzer {
        public:

            HcalTrackAnalyzer(const std::string& name, ldmx::Process& process) : ldmx::Analyzer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void analyze(const ldmx::Event& event);

            virtual void onProcessStart(); 

            virtual void onProcessEnd();

        private:
            
            /** Name of the collection of HcalMipTracks */
            std::string hcalMipTracksCollName_;

            /** Names of the pass that created the HcalMipTracks */
            std::string hcalMipTracksPassName_;
            
            /** Minimum number of non-noise HcalHits for a real track to be findable */
            int minHcalHits_;

            /** Counters for track */
            unsigned int numTracks_[4][4];
           
            /** Number of tracks found per event */
            TH1F* hTracksPerEvent_;

            /** Number of clusters per track */
            TH1F* hClustersPerTrack_;

    };
}

#endif /* HCAL_HCALTRACKANALYZER_H */

