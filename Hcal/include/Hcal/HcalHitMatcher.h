/**
 * @file HcalHitMatcher.h
 * @brief
 * @author Matthew Forsman
 * @author Tom Eichlersmith
 */
#ifndef HCAL_HCALHITMATCHER_H
#define HCAL_HCALHITMATCHER_H

#include "Event/SimTrackerHit.h"
#include "Event/EcalCluster.h"
#include "Event/SimParticle.h"
#include "Framework/EventProcessor.h"
#include "Framework/ParameterSet.h"
#include "Event/CalorimeterHit.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

#include "TH1.h"
#include "TH2.h"
#include "TClonesArray.h"
#include "TVector3.h"

//LDMX Framework
#include "Event/Event.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "Framework/ParameterSet.h" // Needed to import parameters from configuration file

class SimTrackerHit;

namespace ldmx {
    
    /**
     * @class HcalHitMatcher
     * @brief 
     */
    class HcalHitMatcher : public ldmx::Analyzer {
        public:

            HcalHitMatcher(const std::string& name, ldmx::Process& process) : ldmx::Analyzer(name, process) {}

            /**
             * Gets strings from ParameterSet class for collection names.
             * Collection strings that are searched for:
             *      scoringPlaneHits_e (collection of ecal hits that crossed ecal scoring plane)
             *      scoringPlaneHits_h (collection of hcal hits that cross hcal scoring plane)
             */
            virtual void configure(const ldmx::ParameterSet& ps);
            
            /**
             * Fills histograms and attempts to match hcal hits with the SimParticle
             * that caused them.
             */
            virtual void analyze(const ldmx::Event& event);
            
            /**
             * Compares two SimTrackerHit based on the momentum of the SimParticles.
             * Returns true if a has a higher momentum magnitude than b.
             */
            static bool compSimsP(const SimTrackerHit* a, const SimTrackerHit* b);
            
            /**
             * Compares two SimTrackerHits based on the SimParticle.
             * If SimParticles are the same, then compare momentum, otherwise sort by reference.
             */
            static bool compSims(const SimTrackerHit* a, const SimTrackerHit* b); 
    
            /**
             * Calculate the distance between the line segment from v to w and the point p.
             */
            double point_line_distance(TVector3 v, TVector3 w, TVector3 p);
    
            virtual void onFileOpen() { }

            virtual void onFileClose() { }
        
            /**
             * Finds histogram directory and initializes all of the histograms.
             */
            virtual void onProcessStart(); 
            
            /**
             * Prints out totals on numbers of hcal hits.
             */
            virtual void onProcessEnd();

        private:
    
            std::string EcalHitColl_; //* Name of Ecal Digis Collection
            std::string HcalHitColl_; //* Name of Hcal Digis Collection
            std::string EcalScoringPlane_; //* Name of Ecal Scoring Plane Hits Collection
            std::string HcalScoringPlane_; //* Name of Hcal Scoring Plane Hits Collection
            
            long int numNonNoiseHits_; //* Number of Non-Noise Hcal Hits
            long int numMatchedHits_; //* Number of Hcal Hits matched to a sim particle
            long int numEvents_; //* Number of events analyzed
        
            /**
             * Same histograms but limiting particles to specific energy regions
             * These energy regions are specific ranges of standard deviations of the 
             * total particle energy histogram.
             * These ranges are approximate.
             *  Index   Range in Standard Deviations
             *  0       (-1,+1)
             *  1       [+1,+2)
             *  2       (-2,-1]
             *  3       [+2,inf)
             *  4       (-3,-2]
             *  5       (-4,-3]
             *  6       (-5,-4]
             *  7       (-6,-5]
             *  8       (-inf,-6]
             *  9       (-inf,inf) - all events
             */

            //Event information (i.e. One Entry per Event)
            TH1D* h_Ecal_SummedEnergy_SD[10];
            TH1D* h_NumParticles_SD[10];
            TH1D* h_EventMaxPE_SD[10];

            //SimParticle
            TH1D* h_Particle_PDGID_All_SD[10]; //All PDG IDs
            TH1D* h_Particle_PDGID_Matched_SD[10]; //Matched PDG IDs
            TH1D* h_Particle_HitDistance_All_SD[10]; //Distance between SimParticles and HcalHits
            TH1D* h_Particle_HitDistance_Matched_SD[10]; //Distance between SimParticles and HcalHits
            TH1D* h_Particle_Energy_All_SD[10]; //All SimParticle energies
            TH1D* h_Particle_Energy_Matched_SD[10]; //Matched SimParticle energies

            //Position of HcalHits
            TH1D* h_HcalHit_Z_SD[10];
            TH2D* h_HcalHit_ZbyR_All_SD[10];
            TH2D* h_HcalHit_ZbyR_Unmatched_SD[10];
            TH2D* h_HcalHit_ZbyR_TimeLess15_SD[10];
            TH2D* h_HcalHit_ZbyR_TimeGreat40_SD[10];
            TH2D* h_HcalHit_ZbyR_Matched_Photon_SD[10];
            TH2D* h_HcalHit_ZbyR_Matched_Electron_SD[10];
            TH2D* h_HcalHit_ZbyR_Matched_Neutron_SD[10];
            TH2D* h_HcalHit_ZbyR_Matched_Other_SD[10];
            TH2D* h_HcalHit_ZbyR_Matched_TdifLess15_SD[10];
            TH2D* h_HcalHit_ZbyR_Matched_TdifGreat40_SD[10];
            
            //PEs of HcalHit
            TH1D* h_HcalHit_PE_All_SD[10];
            TH1D* h_HcalHit_PE_TimeLess15_SD[10];
            TH1D* h_HcalHit_PE_TimeGreat40_SD[10];
            TH1D* h_HcalHit_PE_Matched_TdifLess15_SD[10];
            TH1D* h_HcalHit_PE_Matched_TdifGreat40_SD[10];

            //Time of HcalHit
            TH1D* h_HcalHit_Time_All_SD[10];
            TH1D* h_HcalHit_Time_Matched_All_SD[10];
            TH1D* h_HcalHit_Time_Matched_Nucleons_SD[10];
            TH1D* h_HcalHit_Time_Matched_Tdif_SD[10];

            //Deprecated for now
//          TH1D* h_HcalHit_photon_energy_SD[10];
//          TH2D* h_HcalHit_nucleon_time_vs_energy_SD[10];

    };
}

#endif /* HCAL_HCALHITMATCHER_H */
