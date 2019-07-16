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

#include "TDirectoryFile.h" //directories to organize histograms

namespace ldmx {

    void HcalHitMatcher::configure(const ldmx::ParameterSet& ps) { 
        
        EcalHitColl_ = ps.getString("EcalHitCollectionName");
        HcalHitColl_ = ps.getString("HcalHitCollectionName");
        EcalScoringPlane_ = ps.getString("EcalScoringPlaneHitsName"); 
        HcalScoringPlane_ = ps.getString("HcalScoringPlaneHitsName");

        return;
    }

    void HcalHitMatcher::analyze(const ldmx::Event& event) {

        //----------This section obtains a list of sim particles that cross the ecal scoring plane---------->
        const TClonesArray* ecalScoringPlaneHits = event.getCollection( EcalScoringPlane_ );
        const TClonesArray* simParticles = event.getCollection("SimParticles"); // side-effect is to make TRefs all valid

        std::vector<SimTrackerHit*> simVec;
        std::vector<SimTrackerHit*> filteredSimVec;//Sim particles are organized from highest to lowest momentum
    
        for (int i = 0; i < ecalScoringPlaneHits->GetEntriesFast(); i++ ) {
            ldmx::SimTrackerHit* ecalSPH = (ldmx::SimTrackerHit*)(ecalScoringPlaneHits->At(i));
            simVec.push_back(ecalSPH);
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

        const TClonesArray* ecalHitColl = event.getCollection( EcalHitColl_ ); 

        double e_cal_sum_energy = 0;
        for(int i=0; i < ecalHitColl->GetEntriesFast(); i++) {
            ldmx::EcalHit* ecalhit = (ldmx::EcalHit*)(ecalHitColl->At(i));
            if ( ! ecalhit->isNoise() ) { //Only add non-noise hits
                e_cal_sum_energy += ecalhit->getEnergy();
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
        const TClonesArray* hcalHitColl = event.getCollection( HcalHitColl_ );

        float max_PE_of_event=0;
        for(int i=0; i < hcalHitColl->GetEntriesFast(); i++) { //Begin loop over hcalhits array
            ldmx::HcalHit* hcalhit = (ldmx::HcalHit*)(hcalHitColl->At(i));
            numHits_++;
            if ( ! hcalhit->getNoise() ) { //Only analyze non-noise hits

                numNonNoiseHits_++;
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
                    TVector3 hCalPoint = TVector3(hcalhit->getX(), hcalhit->getY(), hcalhit->getZ());
                    
                    new_dist = point_line_distance(simStartT, simEndT, hCalPoint);
                    
                    h_ParticleHit_Distance_SD[9]->Fill(new_dist);
                    h_ParticleHit_Distance_SD[ecal_sumESD]->Fill(new_dist);
        
                    if(simStart[2]<10.0 and sP->getEnergy()>3000.0) {
                        //discarding original electron?
                        //check if sim particle is too close to the front of the hcal?
                        //DO NOTHING (skip this sim particle)
                    } else if(new_dist < dist) {
                        dist = new_dist; //Distance to matched particle
                        simPartNum = j; //Matched particle number in array of sim particles
                    }
                    
                } //iterate over sim particles to match one to current hcal hit

                //check if able to match sim particle(s) to hcal hit
                if(simPartNum >= 0) {
                    pdgID = filteredSimVec[simPartNum]->getSimParticle()->getPdgID();
                }

                double hcalhit_radialdist2 = pow(hcalhit->getX(), 2) + pow(hcalhit->getY(), 2);
                double hcalhit_radialdist = 0;
                //check to avoid a floating point error
                if(abs(hcalhit_radialdist2) > 1e-5) {
                    hcalhit_radialdist = sqrt(hcalhit_radialdist2);
                }

                h_HCalhit_zbyr_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                h_ZdepthofHCalHit_SD[9]->Fill(hcalhit->getZ());
                h_hcal_hit_time_all_SD[9]->Fill(hcalhit->getTime());
                h_hcal_hits_all_PEs_SD[9]->Fill(hcalhit->getPE());
                
                h_HCalhit_zbyr_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                h_ZdepthofHCalHit_SD[ecal_sumESD]->Fill(hcalhit->getZ());
                h_hcal_hit_time_all_SD[ecal_sumESD]->Fill(hcalhit->getTime());
                h_hcal_hits_all_PEs_SD[ecal_sumESD]->Fill(hcalhit->getPE());
                
                if(hcalhit->getTime() < 15.0)  {
                    h_hCalhit_time_less15_PE_SD[9]->Fill(hcalhit->getPE());
                    h_hCalhit_time_less15_PE_SD[ecal_sumESD]->Fill(hcalhit->getPE());
                    h_hCalhit_time_less15_position_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                    h_hCalhit_time_less15_position_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                } else if(hcalhit->getTime() > 40.0)  {
                    h_hCalhit_time_great40_PE_SD[9]->Fill(hcalhit->getPE());
                    h_hCalhit_time_great40_PE_SD[ecal_sumESD]->Fill(hcalhit->getPE());
                    h_hCalhit_time_great40_position_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                    h_hCalhit_time_great40_position_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                }
                
                if(hcalhit->getPE() > max_PE_of_event)
                    max_PE_of_event=hcalhit->getPE();
                
                if( dist <= 150.0 ) {//must be 150mm or closer to confidentally match an HCal hit to a sim particle
                
                    numMatchedHits_++;
                    h_HCalhit_getTime_SD[9]->Fill(hcalhit->getTime());
                    h_HCalhit_getTime_SD[ecal_sumESD]->Fill(hcalhit->getTime());
                    
                    double part_hcalhit_timeDiff = (hcalhit->getTime()) - (filteredSimVec[simPartNum]->getSimParticle()->getTime());
                    
                    h_hit_time_creation_time_diff_SD[9]->Fill(part_hcalhit_timeDiff);
                    h_hit_time_creation_time_diff_SD[ecal_sumESD]->Fill(part_hcalhit_timeDiff);
                    if(part_hcalhit_timeDiff < 15.0)  {
                        h_part_hCalhit_tdif_less15_PE_SD[9]->Fill(hcalhit->getPE());
                        h_part_hCalhit_tdif_less15_PE_SD[ecal_sumESD]->Fill(hcalhit->getPE());
                        h_part_hCalhit_tdif_less15_position_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                        h_part_hCalhit_tdif_less15_position_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                    } else if(part_hcalhit_timeDiff > 40.0)  {
                        h_part_hCalhit_tdif_great40_PE_SD[9]->Fill(hcalhit->getPE());
                        h_part_hCalhit_tdif_great40_PE_SD[ecal_sumESD]->Fill(hcalhit->getPE());
                        h_part_hCalhit_tdif_great40_position_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                        h_part_hCalhit_tdif_great40_position_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                    }
                    
                    //protons or neutrons (nucleons) 
                    if( pdgID==2112 or pdgID==2212 ) { 
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
                        case 11:
                            h_HCalhit_electron_zbyr_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            h_HCalhit_electron_zbyr_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            break;
                        case 22:
                            h_HCalhit_photon_zbyr_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            h_HCalhit_photon_zbyr_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                            h_HCalhit_photon_energy_SD[9]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                            h_HCalhit_photon_energy_SD[ecal_sumESD]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                            break;
                        case 2112:
                            h_HCalhit_neutron_zbyr_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                            h_HCalhit_neutron_zbyr_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            break;
                        default:
                            h_HCalhit_other_zbyr_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            h_HCalhit_other_zbyr_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            break;
                    }

                } else {

                    h_HCalhit_unmatched_zbyr_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                    h_HCalhit_unmatched_zbyr_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);

                } //matched or unmatched
    
            } // if not a noise hit

        }//End loop over hcalhits array

        // maximum PE in hcal hits for the event
        h_hcal_hits_max_PE_of_event_SD[9]->Fill(max_PE_of_event);
        h_hcal_hits_max_PE_of_event_SD[ecal_sumESD]->Fill(max_PE_of_event);

        return;
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
    
        //Define a parameter t clamped between 0 and 1 so that the distance is never measured from 
        //beyond the endpoints of the line segment
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
        
        numHits_ = 0;
        numNonNoiseHits_ = 0;
        numMatchedHits_ = 0;

        TDirectory* d_base = getHistoDirectory();

        //create a directory for each type of histogram we care about
        TDirectory* d_PDGIDs = d_base->mkdir( "PDG_IDs" , "PDG_IDs");
        TDirectory* d_ZdeptofHcalHit = d_base->mkdir( "ZDepth_of_HcalHit" , "ZDepth_of_HcalHit" );
        TDirectory* d_ParticleHit_Distance
        TDirectory* d_HCalhit_zbyr
        TDirectory* d_HCalhit_photon_zbyr
        TDirectory* d_HCalhit_electron_zbyr
        TDirectory* d_HCalhit_neutron_zbyr
        TDirectory* d_HCalhit_other_zbyr
        TDirectory* d_HCalhit_unmatched_zbyr
        TDirectory* d_HCalhit_photon_energy
        TDirectory* d_HCalhit_getTime
        TDirectory* d_HCalhit_getTime_nucleons
        TDirectory* d_HCalhit_nucleon_time_vs_energy
        TDirectory* d_E_cal_summed_energy
        TDirectory* d_total_particles
        TDirectory* d_particle_energy
        TDirectory* d_hcal_hit_time_all
        TDirectory* d_hit_time_creation_time_diff
        TDirectory* d_part_hCalhit_tdif_less15_PE
        TDirectory* d_part_hCalhit_tdif_less15_position
        TDirectory* d_part_hCalhit_tdif_great40_PE
        TDirectory* d_part_hCalhit_tdif_great40_position
        TDirectory* d_hCalhit_time_less15_PE
        TDirectory* d_hCalhit_time_less15_position
        TDirectory* d_hCalhit_time_great40_PE
        TDirectory* d_hCalhit_time_great40_position
        TDirectory* d_hcal_hits_all_PEs
        TDirectory* d_hcal_hits_max_PE_of_event

        std::vector<std::string> range_names(10);
        // n <--> negative
        // p <--> positive
        range_names[0] = "n1-p1";
        range_names[1] = "p1-p2";
        range_names[2] = "n2-n1";
        range_names[3] = "p2-pinf";
        range_names[4] = "n3-n2";
        range_names[5] = "n4-n3";
        range_names[6] = "n5-n4";
        range_names[7] = "n6-n5";
        range_names[8] = "ninf-n6]";
        range_names[9] = "ninf-pinf";
   

        for ( int i = 0; i < range_names.size(); i++ ) {
            std::string range_name = range_names.at(i);

            d_PDGIDs->cd();
            h_PDGIDs_SD[i]=new TH1F(
                    ("PDG_IDs_SD_"+range_name).c_str(),
                    ("PDG_IDs_SD_"+range_name).c_str(),
                    10000,-5000,5000);

            d_ZdepthofHcalHit->cd();
            h_ZdepthofHCalHit_SD[i]=new TH1F(
                    ("Z depth of HCal hits_SD_"+range_name).c_str(),
                    ("Z depth of HCal hits_SD"+range_name+" (10mm bins)").c_str(),
                    320, 0, 3200);

            h_ParticleHit_Distance_SD[i]=new TH1F(
                    ("Distance between sim particle and HCal hits_SD_"+range_name).c_str(), 
                    ("Distance between sim particle and HCal hits_SD_"+range_name+" (5mm bins)").c_str(), 
                    400, 0, 2000);
            
            h_HCalhit_zbyr_SD[i]=new TH2D(
                    ("HCal hit locations_SD_"+range_name).c_str(), 
                    ("HCal hit locations_SD_"+range_name+";Z depth (mm); radial distance from z-axis (mm)").c_str(),
                    80,0,3200,112,0,4500);
            h_HCalhit_zbyr_SD[i]->SetDrawOption("colz");

            h_HCalhit_photon_zbyr_SD[i]=new TH2D(
                    ("HCal photon hit locations_SD_"+range_name).c_str(), 
                    ("HCal photon hit locations_SD_"+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);
            h_HCalhit_photon_zbyr_SD[i]->SetDrawOption("colz");

            h_HCalhit_electron_zbyr_SD[i]=new TH2D(
                    ("HCal electron hit locations_SD_"+range_name).c_str(), 
                    ("HCal electron hit locations_SD_"+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);
            h_HCalhit_electron_zbyr_SD[i]->SetDrawOption("colz");

            h_HCalhit_neutron_zbyr_SD[i]=new TH2D(
                    ("HCal neutron hit locations_SD_"+range_name).c_str(), 
                    ("HCal neutron hit locations_SD"+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);
            h_HCalhit_neutron_zbyr_SD[i]->SetDrawOption("colz");

            h_HCalhit_other_zbyr_SD[i]=new TH2D(
                    ("HCal other particle hit locations_SD_"+range_name).c_str(), 
                    ("HCal other particle hit locations_SD_"+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);
            h_HCalhit_other_zbyr_SD[i]->SetDrawOption("colz");

            h_HCalhit_unmatched_zbyr_SD[i]=new TH2D(
                    ("HCal unmatched hit locations_SD_"+range_name).c_str(), 
                    ("HCal unmatched hit locations_SD_"+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);
            h_HCalhit_unmatched_zbyr_SD[i]->SetDrawOption("colz");

            h_HCalhit_photon_energy_SD[i]=new TH1F(
                    ("HCal Photon hit Energies_SD_"+range_name).c_str(),
                    ("HCal Photon hit Energies_SD_"+range_name+";Energy(MeV);Count").c_str(), 
                    4000, 0, 4000);

            h_HCalhit_getTime_SD[i]=new TH1F(
                    ("Creation time of particles causing HCal hits_SD_"+range_name).c_str(), 
                    ("Creation time of particles causing HCal hits_SD_"+range_name+";Time(ns);Number of particles created").c_str(), 
                    500, 0, 500);

            h_HCalhit_getTime_nucleons_SD[i]=new TH1F(
                    ("Creation time of nucleons causing HCal hits_SD_"+range_name).c_str(),
                    ("nucleons causing HCal hits_SD_"+range_name+";Time(ns);Number of Nucleons created").c_str(), 
                    500, 0, 500);

            h_HCalhit_nucleon_time_vs_energy_SD[i]=new TH2D(
                    ("Nucleon time vs energy_SD_"+range_name).c_str(), 
                    ("Nucleon time vs energy_SD_"+range_name+";Creation Time(ns);Energy(MeV)").c_str(),
                    100,0,500,250,0,4000); // 5ns time resolution, 16MeV resolution
            h_HCalhit_nucleon_time_vs_energy_SD[i]->SetDrawOption("colz");

            h_E_cal_summed_energy_SD[i]=new TH1F(
                    ("E_cal_summed_energy_SD_"+range_name).c_str(),
                    ("E_cal_summed_energy_SD_"+range_name+";Energy(MeV)(10MeV bin width);Count").c_str(),
                    800,0,8000);//10MeV bins
            h_E_cal_summed_energy_SD[i]->SetDrawOption("colz");

            h_total_particles_SD[i]=new TH1F(
                    ("total_particles_SD_"+range_name).c_str(),
                    ("total_particles_SD_"+range_name+";Number of particles per event;Count").c_str(),
                    50,0,50);

            h_particle_energy_SD[i]=new TH1F(
                    ("matched_particle_energy_SD_"+range_name).c_str(),
                    ("matched_particle_energy_SD_"+range_name+";Energy(MeV)(5MeV bin width);Count").c_str(),
                    800,0,4000);

            h_hcal_hit_time_all_SD[i]=new TH1F(
                    ("HCal_hit_time_all_SD_"+range_name).c_str(),
                    ("HCal_hit_time_all_SD_"+range_name+";time(ns)(5ns bin width);Count").c_str(),
                    100,0,500);

            h_hit_time_creation_time_diff_SD[i]=new TH1F(
                    ("hit_time_creation_time_diff_SD_"+range_name).c_str(),
                    ("hit_time_creation_time_diff_SD_"+range_name+";time(ns)(2ns bin width);Count").c_str(),
                    100,0,200);

            h_part_hCalhit_tdif_less15_PE_SD[i]=new TH1F(
                    ("part_hCalhit_tdif_less15_PE_SD_"+range_name).c_str(),
                    ("part_hCalhit_tdif_less15_PE_SD_"+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            h_part_hCalhit_tdif_less15_position_SD[i]=new TH2D(
                    ("part_hCalhit_tdif_less15_position_SD_"+range_name).c_str(),
                    ("part_hCalhit_tdif_less15_position_SD_"+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);
            h_part_hCalhit_tdif_less15_position_SD[i]->SetDrawOption("colz");

            h_part_hCalhit_tdif_great40_PE_SD[i]=new TH1F(
                    ("part_hCalhit_tdif_great40_PE_SD_"+range_name).c_str(),
                    ("part_hCalhit_tdif_great40_PE_SD_"+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            h_part_hCalhit_tdif_great40_position_SD[i]=new TH2D(
                    ("part_hCalhit_tdif_great40_position_SD_"+range_name).c_str(),
                    ("part_hCalhit_tdif_great40_position_SD_"+range_name).c_str(),
                    80,0,3200,112,0,4500);
            h_part_hCalhit_tdif_great40_position_SD[i]->SetDrawOption("colz");

            h_hCalhit_time_less15_PE_SD[i]=new TH1F(
                    ("hCalhit_time_less15_PE_SD_"+range_name).c_str(),
                    ("hCalhit_time_less15_PE_SD_"+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            h_hCalhit_time_less15_position_SD[i]=new TH2D(
                    ("hCalhit_time_less15_position_SD_"+range_name).c_str(),
                    ("hCalhit_time_less15_position_SD_"+range_name).c_str(),
                    80,0,3200,112,0,4500);
            h_hCalhit_time_less15_position_SD[i]->SetDrawOption("colz");

            h_hCalhit_time_great40_PE_SD[i]=new TH1F(
                    ("hCalhit_time_great40_PE_SD_"+range_name).c_str(),
                    ("hCalhit_time_great40_PE_SD_"+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            h_hCalhit_time_great40_position_SD[i]=new TH2D(
                    ("hCalhit_time_great40_position_SD_"+range_name).c_str(),
                    ("hCalhit_time_great40_position_SD_"+range_name).c_str(),
                    80,0,3200,112,0,4500);
            h_hCalhit_time_great40_position_SD[i]->SetDrawOption("colz");

            h_hcal_hits_all_PEs_SD[i]=new TH1F(
                    ("hcal_hits_all_PEs_SD_"+range_name).c_str(),
                    ("hcal_hits_all_PEs_SD_"+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            h_hcal_hits_max_PE_of_event_SD[i]=new TH1F(
                    ("h_hcal_hits_max_PE_of_event_SD_"+range_name).c_str(),
                    ("h_hcal_hits_max_PE_of_event_SD_"+range_name+";Photoelectrons(PEs);Count").c_str(),
                    500,0,500);
        } //iterate over range_names

        //go back to base histogram directory
        d_base->cd();
    
        return;
    } //onProcessStart
    
    void HcalHitMatcher::onProcessEnd() {
        
        std::cout << "Total Number of Hits:     " << numHits_ << std::endl;
        std::cout << "Number of Non Noise Hits: " << numNonNoiseHits_ << std::endl;
        std::cout << "Number of Matched Hits:   " << numMatchedHits_ << std::endl;

        return;
    }

} //ldmx namespace

DECLARE_ANALYZER_NS(ldmx, HcalHitMatcher);
