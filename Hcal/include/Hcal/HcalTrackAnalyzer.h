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
#include <set> //std::set of hitkeys from tracks

//ROOT
#include "TH1.h" //One Dimensional Histograms
#include "TClonesArray.h" //Array to get list of hits from Event

//LDMX Framework
#include "Framework/EventProcessor.h" //Needed to declare analyzer
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file
#include "Event/Event.h" //Study event by event
#include "Event/HcalTrack.h" //Study hcal tracks

namespace ldmx {
    
    /**
     * @class HcalTrackAnalyzer
     * @brief ldmx::Analyzer that constructs histograms studying how tracks in the Hcal behave differently.
     *
     * @note Currently, the Hcal strip orientation is not specified (x or y), and simulation is only split
     *  along y (as of June 18, 2018). This means any orientation-related studies in this analyzer will be
     *  making external assumptions that SHOULD BE REMOVED if strip orientation is specifed in the future.
     */
    class HcalTrackAnalyzer : public ldmx::Analyzer {
        public:

            HcalTrackAnalyzer(const std::string& name, ldmx::Process& process) : ldmx::Analyzer(name, process) {}

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void analyze(const ldmx::Event& event);

            virtual void onFileOpen() {}

            virtual void onFileClose() {}

            virtual void onProcessStart(); 

            virtual void onProcessEnd() {}

        private:
            
            void checkForNullPtr( HcalTrack *track ) const;
            
            std::string trackcollname_; //* name of the collection of HcalTracks 
            
            TH1F* h_tracksperevent_; //* number of tracks per event
            TH1F* h_layhitspertrack_[3]; //* number of layer hits per track, index i is the ith track found in the event
    };
}

#endif /* HCAL_HCALTRACKANALYZER_H */

