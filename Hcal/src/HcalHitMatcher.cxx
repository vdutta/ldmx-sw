/**
 * @file HcalHitMatcher.cxx
 * @brief The purpose of this analyzer is to study vetoes caused by activity in the Hcal using 
 *        Monte Carlo simulations.
 *        It extracts an array of sim particles and then matches each sim particle to a hit in the Hcal.
 *        Hcal hits and sim particles are matched spatially by finding the closest distance from a 
 *        sim particle's trajectory to a reconstructed Hcal hit.
 *        Plots and results are then tabulated in ROOT based on the sim particle and Hcal hit matches.
 *
 * @author Matthew Forsman 
 */

#include "Hcal/HcalHitMatcher.h"
#include "Event/HcalHit.h"
#include "Event/SimParticle.h"

namespace ldmx {

    void HcalHitMatcher::configure(const ldmx::ParameterSet& ps) { 
        caloCol_e= ps.getString("caloHitCollection_e");
        caloCol_h= ps.getString("caloHitCollection_h");
        scoringPlane_e = ps.getString("scoringPlaneHits_e"); 
        scoringPlane_h = ps.getString("scoringPlaneHits_h");

        return;
    }

    void HcalHitMatcher::analyze(const ldmx::Event& event) {

        //----------This section obtains a list of sim particles that cross the ecal scoring plane---------->
        const TClonesArray* scoringPlaneHits=event.getCollection(scoringPlane_e);
        const TClonesArray* simParticles=event.getCollection("SimParticles"); // side-effect is to make TRefs all valid
        
        ldmx::SimTrackerHit* ecalSPP;
        std::vector<SimTrackerHit*> simVec;
        std::vector<SimTrackerHit*> filteredSimVec;//Sim particles are organized from highest to lowest momentum
    
        if (scoringPlaneHits->GetEntriesFast() > 0 ) {
            for (TIter next(scoringPlaneHits); ecalSPP = (ldmx::SimTrackerHit*)next();) {
                simVec.push_back(ecalSPP);
            }
        }
    
        std::sort(simVec.begin(), simVec.end(), compSims);
    
        SimParticle* lastP = 0; //sometimes multiple SP hits from same particle
        for (int j = 0; j < simVec.size(); j++) {
            SimParticle* sP = simVec[j]->getSimParticle();
            if (sP == lastP) continue;
            lastP = sP;
            filteredSimVec.push_back(simVec[j]);
        }
    
        std::sort(filteredSimVec.begin(), filteredSimVec.end(), compSimsP);
    
        //----------This section obtains a list of sim particles that cross the hcal scoring plane---------->
        //Currently it is NOT being used to obtain any information
//        
//        const TClonesArray* scoringPlaneHits_h=event.getCollection(scoringPlane_h);
//        const TClonesArray* simParticles_h=event.getCollection("SimParticles"); // side-effect is to make TRefs all valid
//    
//        ldmx::SimTrackerHit* hcalSPP;
//        std::vector<SimTrackerHit*> simVec_h;
//        std::vector<SimTrackerHit*> filteredSimVec_h;//Sim particles are organized from highest to lowest momentum
//    
//        if ( scoringPlaneHits_h->GetEntriesFast() > 0 ) {
//            for (TIter next(scoringPlaneHits_h); hcalSPP = (ldmx::SimTrackerHit*)next();) {
//                simVec_h.push_back(hcalSPP);
//            }
//        }
//    
//        std::sort(simVec_h.begin(), simVec_h.end(), compSims);
//
//        SimParticle* lastP_h = 0; // sometimes multiple SP hits from same particle
//        for (int j = 0; j < simVec_h.size(); j++) {
//            SimParticle* sP_h = simVec_h[j]->getSimParticle();
//            if (sP_h == lastP_h) continue;
//            lastP = sP_h;
//            filteredSimVec_h.push_back(simVec_h[j]);
//        }
//        
//        std::sort(filteredSimVec_h.begin(), filteredSimVec_h.end(), compSimsP);
    
    
        //----------This section calculates the energy in the ECal---------->
        //Then uses this energy to set standard deviation range
        const TClonesArray* tca_e=event.getCollection("ecalDigis","recon");
    
        double e_cal_sum_energy = 0;
        for(int i=0; i < tca_e->GetEntriesFast(); i++) {
            ldmx::EcalHit* eCalhit = (ldmx::EcalHit*)(tca_e->At(i));
            if ( ! eCalhit->isNoise() ) { //Only add non-noise hits
                e_cal_sum_energy += eCalhit->getEnergy();
            }
        }
        
        //Classify this event as one of the standard deviation regions
        int ecal_sumESD = 1; //Labels what approximate standard deviation range the summed ecal energy is
        if(e_cal_sum_energy<4400 && e_cal_sum_energy>3600)       ecal_sumESD=0;//Normal range (within ~1SD)
        else if(e_cal_sum_energy>=4400 && e_cal_sum_energy<4800) ecal_sumESD=1;//High Range (within than +1SD to +2SD)
        else if(e_cal_sum_energy<=3600 && e_cal_sum_energy>3200) ecal_sumESD=2;//low range (within -1SD to -2SD)
        else if(e_cal_sum_energy>=4800)                          ecal_sumESD=3;//Higher Range (higher than 2SD)
        else if(e_cal_sum_energy<=3200 && e_cal_sum_energy>2800) ecal_sumESD=4;//lower Range (within -2SD to -3SD)
        else if(e_cal_sum_energy<=2800 && e_cal_sum_energy>2400) ecal_sumESD=5;//very low range (within -3SD to -4SD)
        else if(e_cal_sum_energy<=2400 && e_cal_sum_energy>2000) ecal_sumESD=6;//extremely low range (within -4SD to -5SD)
        else if(e_cal_sum_energy<=2000 && e_cal_sum_energy>1600) ecal_sumESD=7;//super-duper low range (within -5SD to -6SD)
        else if(e_cal_sum_energy<=1600)                          ecal_sumESD=8;//mega-ultra-super-duper low range (Less than -6SD)
        else ecal_sumESD = 1;//shouldn't ever get here
        
        //Bin event information
        h_E_cal_summed_energy_SD[9]->Fill(e_cal_sum_energy);
        h_E_cal_summed_energy_SD[ecal_sumESD]->Fill(e_cal_sum_energy);
    
        h_total_particles_SD[9]->Fill(filteredSimVec.size());
        h_total_particles_SD[ecal_sumESD]->Fill(filteredSimVec.size());
        
        //----This section matches HCal hits to sim particles and records results----->
        const TClonesArray* tca_h = event.getCollection( "hcalDigis" , "recon" );

        float max_PE_of_event=0;
        for(int i=0; i < tca_h->GetEntriesFast(); i++) { //Begin loop over hCalhits array
            ldmx::HcalHit* hCalhit = (ldmx::HcalHit*)(tca_h->At(i));
            
            if ( ! hCalhit->getNoise() ) { //Only analyze non-noise hits
    
                int pdgID=0; //PDG of sim particle that matched this hit
                int simPartNum = -1; //index of sim particle that matched this hcal hit

                double new_dist=9999, dist=9998; //initializing distance variables
                for(int j=0; j<filteredSimVec.size(); j++) { //Iterate over all sim particles and match one to an Hcal hit
                    SimParticle* sP = filteredSimVec[j]->getSimParticle();
                    pdgID = sP->getPdgID();
                    
                    std::vector<double> simStart = sP->getVertex();
                    std::vector<double> simEnd = sP->getEndPoint();
        
                    TVector3 simStartT = TVector3(simStart[0], simStart[1], simStart[2]);
                    TVector3 simEndT = TVector3(simEnd[0], simEnd[1], simEnd[2]);
                    TVector3 hCalPoint = TVector3(hCalhit->getX(), hCalhit->getY(), hCalhit->getZ());
                    
                    new_dist = point_line_distance(simStartT, simEndT, hCalPoint);
                    
                    h_ParticleHit_Distance_SD[9]->Fill(new_dist);
                    h_ParticleHit_Distance_SD[ecal_sumESD]->Fill(new_dist);
        
                    if(simStart[2]<10.0 && sP->getEnergy()>3000.0);
                    else if(new_dist < dist) {
                        dist = new_dist; //Distance to matched particle
                        simPartNum = j; //Matched particle number in array of sim particles
                    }
                    
                } //iterate over sim particles to match one to current hcal hit

                //check if able to match sim particle(s) to hcal hit
                if(simPartNum < 0) {
                    pdgID = filteredSimVec[simPartNum]->getSimParticle()->getPdgID();
                }
        
                double hCalhit_radialdist2 = pow(hCalhit->getX(), 2) + pow(hCalhit->getY(), 2);
                double hCalhit_radialdist = 0;
                if(abs(hCalhit_radialdist2) > 1e-5) //check to avoid a floating point error
                    hCalhit_radialdist = sqrt(hCalhit_radialdist2);
                
                h_HCalhit_zbyr_SD[9]->Fill(hCalhit->getZ(), hCalhit_radialdist);
                h_ZdepthofHCalHit_SD[9]->Fill(hCalhit->getZ());
                h_hcal_hit_time_all_SD[9]->Fill(hCalhit->getTime());
                h_hcal_hits_all_PEs_SD[9]->Fill(hCalhit->getPE());
                
                h_HCalhit_zbyr_SD[ecal_sumESD]->Fill(hCalhit->getZ(), hCalhit_radialdist);
                h_ZdepthofHCalHit_SD[ecal_sumESD]->Fill(hCalhit->getZ());
                h_hcal_hit_time_all_SD[ecal_sumESD]->Fill(hCalhit->getTime());
                h_hcal_hits_all_PEs_SD[ecal_sumESD]->Fill(hCalhit->getPE());
                
                if(hCalhit->getTime() < 15.0)  {
                    h_hCalhit_time_less15_PE_SD[9]->Fill(hCalhit->getPE());
                    h_hCalhit_time_less15_PE_SD[ecal_sumESD]->Fill(hCalhit->getPE());
                    h_hCalhit_time_less15_position_SD[9]->Fill(hCalhit->getZ(), hCalhit_radialdist);
                    h_hCalhit_time_less15_position_SD[ecal_sumESD]->Fill(hCalhit->getZ(), hCalhit_radialdist);
                } else if(hCalhit->getTime() > 40.0)  {
                    h_hCalhit_time_great40_PE_SD[9]->Fill(hCalhit->getPE());
                    h_hCalhit_time_great40_PE_SD[ecal_sumESD]->Fill(hCalhit->getPE());
                    h_hCalhit_time_great40_position_SD[9]->Fill(hCalhit->getZ(), hCalhit_radialdist);
                    h_hCalhit_time_great40_position_SD[ecal_sumESD]->Fill(hCalhit->getZ(), hCalhit_radialdist);
                }
                
                if(hCalhit->getPE() > max_PE_of_event)
                    max_PE_of_event=hCalhit->getPE();
        
                if( dist <= 150.0 ) {//must be 150mm or closer to confidentally match an HCal hit to a sim particle
                    h_HCalhit_getTime_SD[9]->Fill(hCalhit->getTime());
                    h_HCalhit_getTime_SD[ecal_sumESD]->Fill(hCalhit->getTime());
                    
                    double part_hCalhit_timeDiff = (hCalhit->getTime()) - (filteredSimVec[simPartNum]->getSimParticle()->getTime());
                    
                    h_hit_time_creation_time_diff->Fill(part_hCalhit_timeDiff);
                    h_hit_time_creation_time_diff_SD[ecal_sumESD]->Fill(part_hCalhit_timeDiff);
                    if(part_hCalhit_timeDiff < 15.0)  {
                        h_part_hCalhit_tdif_less15_PE_SD[9]->Fill(hCalhit->getPE());
                        h_part_hCalhit_tdif_less15_PE_SD[ecal_sumESD]->Fill(hCalhit->getPE());
                        h_part_hCalhit_tdif_less15_position_SD[9]->Fill(hCalhit->getZ(), hCalhit_radialdist);
                        h_part_hCalhit_tdif_less15_position_SD[ecal_sumESD]->Fill(hCalhit->getZ(), hCalhit_radialdist);
                    } else if(part_hCalhit_timeDiff > 40.0)  {
                        h_part_hCalhit_tdif_great40_PE_SD[9]->Fill(hCalhit->getPE());
                        h_part_hCalhit_tdif_great40_PE_SD[ecal_sumESD]->Fill(hCalhit->getPE());
                                        h_part_hCalhit_tdif_great40_position->Fill(hCalhit->getZ(), hCalhit_radialdist);
                        h_part_hCalhit_tdif_great40_position_SD[ecal_sumESD]->Fill(hCalhit->getZ(), hCalhit_radialdist);
                    }
        
                    if(pdgID==2112 || pdgID==2212) { 
                        h_HCalhit_getTime_nucleons_SD[9]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getTime());
                        h_HCalhit_nucleon_time_vs_energy_SD[9]->Fill(
                                filteredSimVec[simPartNum]->getSimParticle()->getTime(),
                                filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                        h_HCalhit_getTime_nucleons_SD[ecal_sumESD]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getTime());
                        h_HCalhit_nucleon_time_vs_energy_SD[ecal_sumESD]->Fill(
                                filteredSimVec[simPartNum]->getSimParticle()->getTime(),
                                filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                    }
                    
                    h_PDGIDs_SD[9]->Fill(pdgID);
                    h_PDGIDs_SD[ecal_sumESD]->Fill(pdgID);
        
                    h_particle_energy_SD[9]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                    h_particle_energy_SD[ecal_sumESD]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
        
                    switch(pdgID) {
                        case 11:h_HCalhit_electron_zbyr_SD[9]->Fill(hCalhit->getZ(), hCalhit_radialdist); 
                            h_HCalhit_electron_zbyr_SD[ecal_sumESD]->Fill(hCalhit->getZ(), hCalhit_radialdist); 
                            break;
                        case 22:h_HCalhit_photon_zbyr_SD[9]->Fill(hCalhit->getZ(), hCalhit_radialdist); 
                            h_HCalhit_photon_zbyr_SD[ecal_sumESD]->Fill(hCalhit->getZ(), hCalhit_radialdist);
                            h_HCalhit_photon_energy_SD[9]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                            h_HCalhit_photon_energy_SD[ecal_sumESD]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                            break;
                        case 2112:h_HCalhit_neutron_zbyr_SD[ecal_sumESD]->Fill(hCalhit->getZ(), hCalhit_radialdist);
                              h_HCalhit_neutron_zbyr_SD[9]->Fill(hCalhit->getZ(), hCalhit_radialdist); 
                              break;
                        default:h_HCalhit_other_zbyr_SD[9]->Fill(hCalhit->getZ(), hCalhit_radialdist); 
                            h_HCalhit_other_zbyr_SD[ecal_sumESD]->Fill(hCalhit->getZ(), hCalhit_radialdist); 
                            break;
                    }
                } else {
                    h_HCalhit_unmatched_zbyr_SD[9]->Fill(hCalhit->getZ(), hCalhit_radialdist);
                    h_HCalhit_unmatched_zbyr_SD[ecal_sumESD]->Fill(hCalhit->getZ(), hCalhit_radialdist);
                } //matched or unmatched
    
            } // if not a noise hit

        }//End loop over hCalhits array

        // maximum PE in hcal hits for the event
        h_hcal_hits_max_PE_of_event_SD[9]->Fill(max_PE_of_event);
        h_hcal_hits_max_PE_of_event_SD[ecal_sumESD]->Fill(max_PE_of_event);

    } //analyze

   
    bool HcalHitMatcher::compSimsP(const SimTrackerHit* a, const SimTrackerHit* b) {
        std::vector<double> paVec = a->getMomentum();
        std::vector<double> pbVec = b->getMomentum();
    
        double pa2 = pow(paVec[0],2)+pow(paVec[1],2)+pow(paVec[2],2);
        double pb2 = pow(pbVec[0],2)+pow(pbVec[1],2)+pow(pbVec[2],2);
    
        return pa2 > pb2;
    }

    bool HcalHitMatcher::compSims(const SimTrackerHit* a, const SimTrackerHit* b) {

        if (a->getSimParticle() == b->getSimParticle()) {
             return compSimsP(a,b);
        } else {
             return a->getSimParticle() < b->getSimParticle();
        }
    }

    //This function finds the minimum distance between a line segment and a point in 3D space
    double HcalHitMatcher::point_line_distance(TVector3 v, TVector3 w, TVector3 p)  {
        
        //Define 2 line segments from point v to w and from point v to p
        TVector3 vw = TVector3(w.x()-v.x(), w.y()-v.y(), w.z()-v.z());
        TVector3 vp = TVector3(v.x()-p.x(), v.y()-p.y(), v.z()-p.z());
    
        //Define both line segment's magnitude squared
        const double L2vw = pow(vw.Mag(), 2);
        const double L2vp = pow(vp.Mag(), 2);
    
        //Define a parameter t clamped between 0 and 1 so that the distance is never measured from beyond the endpoints of the line segment
        const double t = std::max(0.0, -1*vp.Dot(vw)/L2vw);//std::max(0.0, std::min(1.0, -1*vp.Dot(vw)/L2vw));
    
        //Calculate the minimum distance squared
        double d2 = L2vp + 2.0*t*vp.Dot(vw) + pow(t,2)*L2vw;
    
        //Return minimum distance and catch floating point errors near zero
        if(abs(d2) < 1e-5)
            return 0;
        else
            return sqrt(d2);
    }

    void HcalHitMatcher::onProcessStart() {
        
        getHistoDirectory();
        
        std::vector<std::string> range_names(10);
        range_names[0] = "(-1,+1)";
        range_names[1] = "[+1,+2)";
        range_names[2] = "(-2,-1)";
        range_names[3] = "[+2,+inf)";
        range_names[4] = "(-3,-2]";
        range_names[5] = "(-4,-3]";
        range_names[6] = "(-5,-4]";
        range_names[7] = "(-6,-5]";
        range_names[8] = "(-inf,-6]";
        range_names[9] = "(-inf,+inf)";
    
        for ( int i = 0; i < range_names.size(); i++ ) {
            std::string range_name = range_names.at(i);
            h_PDGIDs_SD[i]=new TH1F("PDG_IDs_SD_"+range_name,
                    "PDG_IDs_SD1",10000,-5000,5000);
            h_ZdepthofHCalHit_SD[i]=new TH1F("Z depth of HCal hits_SD_"+range_name, 
                    "Z depth of HCal hits_SD"+range_name+" (10mm bins)", 320, 0, 3200);
            h_ParticleHit_Distance_SD[0]=new TH1F("Distance between sim particle and HCal hits_SD_"+range_name, 
                    "Distance between sim particle and HCal hits_SD_"+range_name+" (5mm bins)", 400, 0, 2000);
            h_HCalhit_zbyr_SD[i]=new TH2D("HCal hit locations_SD_"+range_name, 
                    "HCal hit locations_SD_"+range_name+";Z depth (mm); radial distance from z-axis (mm)",
                    80,0,3200,112,0,4500);
            h_HCalhit_photon_zbyr_SD[i]=new TH2D("HCal photon hit locations_SD_"+range_name, 
                    "HCal photon hit locations_SD_"+range_name+";Z depth(mm);radial distance from z-axis(mm)",
                    80,0,3200,112,0,4500);
            h_HCalhit_electron_zbyr_SD[i]=new TH2D("HCal electron hit locations_SD_"+range_name, 
                    "HCal electron hit locations_SD_"+range_name+";Z depth(mm);radial distance from z-axis(mm)",
                    80,0,3200,112,0,4500);
            h_HCalhit_neutron_zbyr_SD[i]=new TH2D("HCal neutron hit locations_SD_"+range_name, 
                    "HCal neutron hit locations_SD"+range_name+";Z depth(mm);radial distance from z-axis(mm)",
                    80,0,3200,112,0,4500);
            h_HCalhit_other_zbyr_SD[i]=new TH2D("HCal other particle hit locations_SD_"+range_name, 
                    "HCal other particle hit locations_SD_"+range_name+";Z depth(mm);radial distance from z-axis(mm)",
                    80,0,3200,112,0,4500);
            h_HCalhit_unmatched_zbyr_SD[i]=new TH2D("HCal unmatched hit locations_SD_"+range_name, 
                    "HCal unmatched hit locations_SD_"+range_name+";Z depth(mm);radial distance from z-axis(mm)",
                    80,0,3200,112,0,4500);
            h_HCalhit_photon_energy_SD[i]=new TH1F("HCal Photon hit Energies_SD_"+range_name,
                    "HCal Photon hit Energies_SD_"+range_name+";Energy(MeV);Count", 4000, 0, 4000);
            h_HCalhit_getTime_SD[i]=new TH1F("Creation time of particles causing HCal hits_SD_"+range_name, 
                    "Creation time of particles causing HCal hits_SD_"+range_name+";Time(ns);Number of particles created", 500, 0, 500);
            h_HCalhit_getTime_nucleons_SD[i]=new TH1F("Creation time of nucleons causing HCal hits_SD_"+range_name,
                    "nucleons causing HCal hits_SD_"+range_name+";Time(ns);Number of Nucleons created", 500, 0, 500);
            h_HCalhit_nucleon_time_vs_energy_SD[i]=new TH2D("Nucleon time vs energy_SD_"+range_name, 
                    "Nucleon time vs energy_SD_"+range_name+";Creation Time(ns);Energy(MeV)",100,0,500,250,0,4000); // 5ns time resolution, 16MeV resolution
            h_E_cal_summed_energy_SD[i]=new TH1F("E_cal_summed_energy_SD_"+range_name,
                    "E_cal_summed_energy_SD_"+range_name+";Energy(MeV)(10MeV bin width);Count",800,0,8000);//10MeV bins
            h_total_particles_SD[i]=new TH1F("total_particles_SD_"+range_name,
                    "total_particles_SD_+"range_name+";Number of particles per event;Count",50,0,50);
            h_particle_energy_SD[i]=new TH1F("matched_particle_energy_SD_"+range_name,
                    "matched_particle_energy_SD_"+range_name+";Energy(MeV)(5MeV bin width);Count",800,0,4000);
            h_hcal_hit_time_all_SD[i]=new TH1F("HCal_hit_time_all_SD_"+range_name,
                    "HCal_hit_time_all_SD_"+range_name+";time(ns)(5ns bin width);Count",100,0,500);
            h_hit_time_creation_time_diff_SD[i]=new TH1F("hit_time_creation_time_diff_SD_"+range_name,
                    "hit_time_creation_time_diff_SD_"+range_name+";time(ns)(2ns bin width);Count",100,0,200);
            h_part_hCalhit_tdif_less15_PE_SD[i]=new TH1F("part_hCalhit_tdif_less15_PE_SD_"+range_name,
                    "part_hCalhit_tdif_less15_PE_SD_"+range_name+";Photoelectrons(PEs);Count",200,0,200);
            h_part_hCalhit_tdif_less15_position_SD[i]=new TH2D("part_hCalhit_tdif_less15_position_SD_"+range_name,
                    "part_hCalhit_tdif_less15_position_SD_"+range_name+";Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
            h_part_hCalhit_tdif_great40_PE_SD[i]=new TH1F("part_hCalhit_tdif_great40_PE_SD_"+range_name,
                    "part_hCalhit_tdif_great40_PE_SD_"+range_name+";Photoelectrons(PEs);Count",200,0,200);
            h_part_hCalhit_tdif_great40_position_SD[i]=new TH2D("part_hCalhit_tdif_great40_position_SD_"+range_name,
                    "part_hCalhit_tdif_great40_position_SD_"+range_name,80,0,3200,112,0,4500);
            h_hCalhit_time_less15_PE_SD[i]=new TH1F("hCalhit_time_less15_PE_SD_"+range_name,
                    "hCalhit_time_less15_PE_SD_"+range_name+";Photoelectrons(PEs);Count",200,0,200);
            h_hCalhit_time_less15_position_SD[i]=new TH2D("hCalhit_time_less15_position_SD_"+range_name,
                    "hCalhit_time_less15_position_SD_"+range_name,80,0,3200,112,0,4500);
            h_hCalhit_time_great40_PE_SD[i]=new TH1F("hCalhit_time_great40_PE_SD_"+range_name,
                    "hCalhit_time_great40_PE_SD_"+range_name+";Photoelectrons(PEs);Count",200,0,200);
            h_hCalhit_time_great40_position_SD[i]=new TH2D("hCalhit_time_great40_position_SD_"+range_name,
                    "hCalhit_time_great40_position_SD_"+range_name,80,0,3200,112,0,4500);
            h_hcal_hits_all_PEs_SD[i]=new TH1F("hcal_hits_all_PEs_SD_"+range_name,
                    "hcal_hits_all_PEs_SD_"+range_name+";Photoelectrons(PEs);Count",200,0,200);
            h_hcal_hits_max_PE_of_event_SD[i]=new TH1F("h_hcal_hits_max_PE_of_event_SD_"+range_name,
                    "h_hcal_hits_max_PE_of_event_SD_"+range_name+";Photoelectrons(PEs);Count",500,0,500);
        } //iterate over range_names
    
        return;
    } //onProcessStart
    
} //ldmx namespace

DECLARE_ANALYZER_NS(ldmx, HcalHitMatcher);
