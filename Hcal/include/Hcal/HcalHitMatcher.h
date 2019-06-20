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

            virtual void configure(const ldmx::ParameterSet& ps);

            virtual void analyze(const ldmx::Event& event);
            
	    static bool compSimsP(const SimTrackerHit* a, const SimTrackerHit* b);
        
	    static bool compSims(const SimTrackerHit* a, const SimTrackerHit* b); 

	    static double point_line_distance(TVector3 v, TVector3 w, TVector3 p);

	    virtual void onFileOpen();

            virtual void onFileClose();
		
            virtual void onProcessStart(); 

            virtual void onProcessEnd();

        private:
	
	std::string caloCol_e, caloCol_h, scoringPlane_e, scoringPlane_h; //scoringPlane_;
	
	//All events (without parsing by ecal energy SD region)
	TH1* h_PDGIDs;
	TH1* h_ZdepthofHCalHit;
	TH1* h_ParticleHit_Distance;
	TH2D* h_HCalhit_zbyr;
	TH2D* h_HCalhit_photon_zbyr;
	TH2D* h_HCalhit_electron_zbyr;
	TH2D* h_HCalhit_neutron_zbyr;
	TH2D* h_HCalhit_other_zbyr;
	TH2D* h_HCalhit_unmatched_zbyr;
	TH1* h_HCalhit_photon_energy;
	TH1* h_HCalhit_getTime;
	TH1* h_HCalhit_getTime_nucleons;
	TH2D* h_HCalhit_nucleon_time_vs_energy;
	TH1F* h_E_cal_summed_energy;
	TH1F* h_total_particles;
	TH1F* h_particle_energy;
	TH1F* h_hcal_hit_time_all;
	TH1F* h_hit_time_creation_time_diff;
	TH1F* h_part_hCalhit_tdif_less15_PE;
	TH2D* h_part_hCalhit_tdif_less15_position;
	TH1F* h_part_hCalhit_tdif_great40_PE;
	TH2D* h_part_hCalhit_tdif_great40_position;
	TH1F* h_hCalhit_time_less15_PE;
	TH2D* h_hCalhit_time_less15_position;
	TH1F* h_hCalhit_time_great40_PE;
	TH2D* h_hCalhit_time_great40_position;
	TH1F* h_hcal_hits_all_PEs;
	TH1F* h_hcal_hits_max_PE_of_event;

	//SD region 1 (normal range, within ~1SD)
	TH1* h_PDGIDs_SD[9];
	TH1* h_ZdepthofHCalHit_SD[9];
	TH1* h_ParticleHit_Distance_SD[9];
	TH2D* h_HCalhit_zbyr_SD[9];
	TH2D* h_HCalhit_photon_zbyr_SD[9];
	TH2D* h_HCalhit_electron_zbyr_SD[9];
	TH2D* h_HCalhit_neutron_zbyr_SD[9];
	TH2D* h_HCalhit_other_zbyr_SD[9];
	TH2D* h_HCalhit_unmatched_zbyr_SD[9];
	TH1* h_HCalhit_photon_energy_SD[9];
	TH1* h_HCalhit_getTime_SD[9];
	TH1* h_HCalhit_getTime_nucleons_SD[9];
	TH2D* h_HCalhit_nucleon_time_vs_energy_SD[9];
	TH1F* h_E_cal_summed_energy_SD[9];
	TH1F* h_total_particles_SD[9];
	TH1F* h_particle_energy_SD[9];
	TH1F* h_hcal_hit_time_all_SD[9];
	TH1F* h_hit_time_creation_time_diff_SD[9];
	TH1F* h_part_hCalhit_tdif_less15_PE_SD[9];
	TH2D* h_part_hCalhit_tdif_less15_position_SD[9];
	TH1F* h_part_hCalhit_tdif_great40_PE_SD[9];
	TH2D* h_part_hCalhit_tdif_great40_position_SD[9];
	TH1F* h_hCalhit_time_less15_PE_SD[9];
	TH2D* h_hCalhit_time_less15_position_SD[9];
	TH1F* h_hCalhit_time_great40_PE_SD[9];
	TH2D* h_hCalhit_time_great40_position_SD[9];
	TH1F* h_hcal_hits_all_PEs_SD[9];
	TH1F* h_hcal_hits_max_PE_of_event_SD[9];

	};
}

#endif /* HCAL_HCALHITMATCHER_H */
