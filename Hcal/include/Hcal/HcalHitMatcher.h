/**
 * @file HcalHitMatcher.h
 * @brief
 * @author Matthew Forsman
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
             *      caloHitCollection_e (collection of ecal hits)
             *      caloHitCollection_h (collection of hcal hits)
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

            virtual void onProcessEnd() { }

        private:
	
        	std::string caloCol_e, caloCol_h, scoringPlane_e, scoringPlane_h; //scoringPlane_;
        
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
        	TH1* h_PDGIDs_SD[10];
        	TH1* h_ZdepthofHCalHit_SD[10];
        	TH1* h_ParticleHit_Distance_SD[10];
        	TH2D* h_HCalhit_zbyr_SD[10];
        	TH2D* h_HCalhit_photon_zbyr_SD[10];
        	TH2D* h_HCalhit_electron_zbyr_SD[10];
        	TH2D* h_HCalhit_neutron_zbyr_SD[10];
        	TH2D* h_HCalhit_other_zbyr_SD[10];
        	TH2D* h_HCalhit_unmatched_zbyr_SD[10];
        	TH1* h_HCalhit_photon_energy_SD[10];
        	TH1* h_HCalhit_getTime_SD[10];
        	TH1* h_HCalhit_getTime_nucleons_SD[10];
        	TH2D* h_HCalhit_nucleon_time_vs_energy_SD[10];
        	TH1F* h_E_cal_summed_energy_SD[10];
        	TH1F* h_total_particles_SD[10];
        	TH1F* h_particle_energy_SD[10];
        	TH1F* h_hcal_hit_time_all_SD[10];
        	TH1F* h_hit_time_creation_time_diff_SD[10];
        	TH1F* h_part_hCalhit_tdif_less15_PE_SD[10];
        	TH2D* h_part_hCalhit_tdif_less15_position_SD[10];
        	TH1F* h_part_hCalhit_tdif_great40_PE_SD[10];
        	TH2D* h_part_hCalhit_tdif_great40_position_SD[10];
        	TH1F* h_hCalhit_time_less15_PE_SD[10];
        	TH2D* h_hCalhit_time_less15_position_SD[10];
        	TH1F* h_hCalhit_time_great40_PE_SD[10];
        	TH2D* h_hCalhit_time_great40_position_SD[10];
        	TH1F* h_hcal_hits_all_PEs_SD[10];
        	TH1F* h_hcal_hits_max_PE_of_event_SD[10];

	};
}

#endif /* HCAL_HCALHITMATCHER_H */
