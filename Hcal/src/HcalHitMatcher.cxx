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
        
        const TClonesArray* scoringPlaneHits_h=event.getCollection(scoringPlane_h);
        const TClonesArray* simParticles_h=event.getCollection("SimParticles"); // side-effect is to make TRefs all valid
    
        ldmx::SimTrackerHit* hcalSPP;
        std::vector<SimTrackerHit*> simVec_h;
        std::vector<SimTrackerHit*> filteredSimVec_h;//Sim particles are organized from highest to lowest momentum
    
        if ( scoringPlaneHits_h->GetEntriesFast() > 0 ) {
            for (TIter next(scoringPlaneHits_h); hcalSPP = (ldmx::SimTrackerHit*)next();) {
                simVec_h.push_back(hcalSPP);
            }
        }
    
        std::sort(simVec_h.begin(), simVec_h.end(), compSims);

        SimParticle* lastP_h = 0; // sometimes multiple SP hits from same particle
        for (int j = 0; j < simVec_h.size(); j++) {
            SimParticle* sP_h = simVec_h[j]->getSimParticle();
            if (sP_h == lastP_h) continue;
            lastP = sP_h;
            filteredSimVec_h.push_back(simVec_h[j]);
        }
        
        std::sort(filteredSimVec_h.begin(), filteredSimVec_h.end(), compSimsP);
    
    
        //----------This section calculates the energy in the ECal---------->
        int ecalHitCount = 0;
        const TClonesArray* tca_e=event.getCollection("ecalDigis","recon");
        const ldmx::EcalHit* eCalhits[tca_e->GetEntriesFast()];
    
        for(size_t i=0; i<tca_e->GetEntriesFast(); i++) {
            ldmx::EcalHit* chit_e = (ldmx::EcalHit*)(tca_e->At(i));
            if(chit_e->isNoise()) continue;
            else {
                eCalhits[ecalHitCount] = chit_e;
                ecalHitCount++;
            }
        }
    
        double e_cal_sum_energy = 0;
        for(int i=0; i<ecalHitCount; i++) {
            const ldmx::EcalHit* eCalhit = (const ldmx::EcalHit*)(eCalhits[i]);
            e_cal_sum_energy += eCalhit->getEnergy();
        }
    
        int ecal_sumESD = 1; //Labels what approximate standard deviation range the summed ecal energy is
        if(e_cal_sum_energy<4400 && e_cal_sum_energy>3600) ecal_sumESD=0; //Normal range (within ~1SD)
        else if(e_cal_sum_energy>=4400 && e_cal_sum_energy<4800) ecal_sumESD=1; //High Range (within than +1SD to +2SD)
        else if(e_cal_sum_energy<=3600 && e_cal_sum_energy>3200) ecal_sumESD=2;//low range (within -1SD to -2SD)
        else if(e_cal_sum_energy>=4800) ecal_sumESD=3; //Higher Range (higher than 2SD)
        else if(e_cal_sum_energy<=3200 && e_cal_sum_energy>2800) ecal_sumESD=4;//lower Range (within -2SD to -3SD)
        else if(e_cal_sum_energy<=2800 && e_cal_sum_energy>2400) ecal_sumESD=5;//very low range (within -3SD to -4SD)
        else if(e_cal_sum_energy<=2400 && e_cal_sum_energy>2000) ecal_sumESD=6;//extremely low range (within -4SD to -5SD)
        else if(e_cal_sum_energy<=2000 && e_cal_sum_energy>1600) ecal_sumESD=7;//super-duper low range (within -5SD to -6SD)
        else if(e_cal_sum_energy<=1600) ecal_sumESD=8;//mega-ultra-super-duper low range (Less than -6SD)
        else ecal_sumESD= 1;//shouldn't ever get here
        
        h_E_cal_summed_energy->Fill(e_cal_sum_energy);
        h_E_cal_summed_energy_SD[ecal_sumESD]->Fill(e_cal_sum_energy);
    
        h_total_particles->Fill(filteredSimVec.size());
        h_total_particles_SD[ecal_sumESD]->Fill(filteredSimVec.size());
    
        //----------This section obtains the HCal hits----------> 
        int hcalHitCount = 0;
        const TClonesArray* tca_h=event.getCollection("hcalDigis","recon");
        const ldmx::HcalHit* hCalhits[tca_h->GetEntriesFast()];
    
        for(size_t i=0; i<tca_h->GetEntriesFast(); i++) {
            ldmx::HcalHit* chit_h = (ldmx::HcalHit*)(tca_h->At(i));
            if((chit_h->getX()==0) && (chit_h->getY()==0) && (chit_h->getZ()==0)) continue; //This filters noise hits
            else if((chit_h->getPE()>=2)) {
                hCalhits[hcalHitCount] = chit_h;
                hcalHitCount++;
            }
        }
        
        //----This section matches HCal hits to sim particles and records results----->
        float max_PE_of_event=0;
        for(int i=0; i<hcalHitCount; i++) { //Begin loop over hCalhits array
            int pdgID=0, new_dist=9999, dist=9998, simPartNum = -1; //Probably a better way to record distance variables
            const ldmx::HcalHit* hCalhit = (const ldmx::HcalHit*)(hCalhits[i]);
            for(int j=0; j<filteredSimVec.size(); j++) { //Iterate over all sim particles and match one to an Hcal hit
                SimParticle* sP = filteredSimVec[j]->getSimParticle();
                pdgID = sP->getPdgID();
                
                std::vector<double> simStart = sP->getVertex();
                std::vector<double> simEnd = sP->getEndPoint();
    
                TVector3 simStartT = TVector3(simStart[0], simStart[1], simStart[2]);
                TVector3 simEndT = TVector3(simEnd[0], simEnd[1], simEnd[2]);
                TVector3 hCalPoint = TVector3(hCalhits[i]->getX(), hCalhits[i]->getY(), hCalhits[i]->getZ());
                
                new_dist = point_line_distance(simStartT, simEndT, hCalPoint);
                
                h_ParticleHit_Distance->Fill(new_dist);
                h_ParticleHit_Distance_SD[ecal_sumESD]->Fill(new_dist);
    
                if(simStart[2]<10.0 && sP->getEnergy()>3000.0);
                else if(new_dist<=dist) {
                    dist = new_dist; //Distance to matched particle
                    simPartNum = j; //Matched particle number in array of sim particles
                }
                
            }
            if(simPartNum != -1) pdgID = filteredSimVec[simPartNum]->getSimParticle()->getPdgID();
    
            double hCalhit_radialdist2 = pow(hCalhits[i]->getX(), 2) + pow(hCalhits[i]->getY(), 2);
            double hCalhit_radialdist = 0;
            if(abs(hCalhit_radialdist2) > 1e-5) hCalhit_radialdist = sqrt(hCalhit_radialdist2);//To avoid FP errors
            
            h_HCalhit_zbyr->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
            h_ZdepthofHCalHit->Fill(hCalhits[i]->getZ());
            h_hcal_hit_time_all->Fill(hCalhits[i]->getTime());
            h_hcal_hits_all_PEs->Fill(hCalhits[i]->getPE());
            
            h_HCalhit_zbyr_SD[ecal_sumESD]->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
            h_ZdepthofHCalHit_SD[ecal_sumESD]->Fill(hCalhits[i]->getZ());
            h_hcal_hit_time_all_SD[ecal_sumESD]->Fill(hCalhits[i]->getTime());
            h_hcal_hits_all_PEs_SD[ecal_sumESD]->Fill(hCalhits[i]->getPE());
            
            if(hCalhits[i]->getTime() < 15.0)  {
                h_hCalhit_time_less15_PE->Fill(hCalhits[i]->getPE());
                h_hCalhit_time_less15_PE_SD[ecal_sumESD]->Fill(hCalhits[i]->getPE());
                h_hCalhit_time_less15_position->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
                h_hCalhit_time_less15_position_SD[ecal_sumESD]->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
            }
            else if(hCalhits[i]->getTime() > 40.0)  {
                h_hCalhit_time_great40_PE->Fill(hCalhits[i]->getPE());
                h_hCalhit_time_great40_PE_SD[ecal_sumESD]->Fill(hCalhits[i]->getPE());
                            h_hCalhit_time_great40_position->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
                            h_hCalhit_time_great40_position_SD[ecal_sumESD]->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
            }
            
            if(i==hcalHitCount-1)  {
                if(hCalhits[i]->getPE() > max_PE_of_event) max_PE_of_event=hCalhits[i]->getPE();
                h_hcal_hits_max_PE_of_event->Fill(max_PE_of_event);
                h_hcal_hits_max_PE_of_event_SD[ecal_sumESD]->Fill(max_PE_of_event);
            } else {
                if(hCalhits[i]->getPE() > max_PE_of_event) max_PE_of_event=hCalhits[i]->getPE();
            }
    
            if(dist<=150.0) {//must be 150mm or closer to confidentally match an HCal hit to a sim particle
                h_HCalhit_getTime->Fill(hCalhits[i]->getTime());
                h_HCalhit_getTime_SD[ecal_sumESD]->Fill(hCalhits[i]->getTime());
                
                double part_hCalhit_timeDiff = (hCalhits[i]->getTime()) - (filteredSimVec[simPartNum]->getSimParticle()->getTime());
                
                h_hit_time_creation_time_diff->Fill(part_hCalhit_timeDiff);
                h_hit_time_creation_time_diff_SD[ecal_sumESD]->Fill(part_hCalhit_timeDiff);
                if(part_hCalhit_timeDiff < 15.0)  {
                    h_part_hCalhit_tdif_less15_PE->Fill(hCalhits[i]->getPE());
                    h_part_hCalhit_tdif_less15_PE_SD[ecal_sumESD]->Fill(hCalhits[i]->getPE());
                    h_part_hCalhit_tdif_less15_position->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
                    h_part_hCalhit_tdif_less15_position_SD[ecal_sumESD]->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
                }
                else if(part_hCalhit_timeDiff > 40.0)  {
                    h_part_hCalhit_tdif_great40_PE->Fill(hCalhits[i]->getPE());
                    h_part_hCalhit_tdif_great40_PE_SD[ecal_sumESD]->Fill(hCalhits[i]->getPE());
                                    h_part_hCalhit_tdif_great40_position->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
                    h_part_hCalhit_tdif_great40_position_SD[ecal_sumESD]->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
                }
    
                if(pdgID==2112 || pdgID==2212) { 
                    h_HCalhit_getTime_nucleons->Fill(filteredSimVec[simPartNum]->getSimParticle()->getTime());
                    h_HCalhit_nucleon_time_vs_energy->Fill(filteredSimVec[simPartNum]->getSimParticle()->getTime(), filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                    h_HCalhit_getTime_nucleons_SD[ecal_sumESD]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getTime());
                                    h_HCalhit_nucleon_time_vs_energy_SD[ecal_sumESD]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getTime(), filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                }
                
                h_PDGIDs->Fill(pdgID);
                h_PDGIDs_SD[ecal_sumESD]->Fill(pdgID);
    
                h_particle_energy->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                h_particle_energy_SD[ecal_sumESD]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
    
                switch(pdgID) {
                    case 11:h_HCalhit_electron_zbyr->Fill(hCalhits[i]->getZ(), hCalhit_radialdist); 
                        h_HCalhit_electron_zbyr_SD[ecal_sumESD]->Fill(hCalhits[i]->getZ(), hCalhit_radialdist); 
                        break;
                    case 22:h_HCalhit_photon_zbyr->Fill(hCalhits[i]->getZ(), hCalhit_radialdist); 
                        h_HCalhit_photon_zbyr_SD[ecal_sumESD]->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
                        h_HCalhit_photon_energy->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                        h_HCalhit_photon_energy_SD[ecal_sumESD]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                        break;
                    case 2112:h_HCalhit_neutron_zbyr_SD[ecal_sumESD]->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
                          h_HCalhit_neutron_zbyr->Fill(hCalhits[i]->getZ(), hCalhit_radialdist); 
                          break;
                    default:h_HCalhit_other_zbyr->Fill(hCalhits[i]->getZ(), hCalhit_radialdist); 
                        h_HCalhit_other_zbyr_SD[ecal_sumESD]->Fill(hCalhits[i]->getZ(), hCalhit_radialdist); 
                        break;
                }
    
            }
            else{
                h_HCalhit_unmatched_zbyr->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
                h_HCalhit_unmatched_zbyr_SD[ecal_sumESD]->Fill(hCalhits[i]->getZ(), hCalhit_radialdist);
            }
        }//End loop over hCalhits array

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
        std::cout << "HcalHitMatcher: Starting processing!" << std::endl;
        getHistoDirectory();
        
        //All events (without parsing by energy SD region)
        h_PDGIDs=new TH1F("PDG_IDs","PDG_IDs",10000,-5000,5000);
        h_ZdepthofHCalHit=new TH1F("Z depth of HCal hits(10mm bins)", "Z depth of HCal hits(10mm bins)", 320, 0, 3200);
        h_ParticleHit_Distance=new TH1F("Distance between sim particle and HCal hits(5mm bins)", "Distance between sim particle and HCal hits(5mm bins)", 400, 0, 2000);
        h_HCalhit_zbyr=new TH2D("HCal hit locations", "HCal hit locations;Z depth (mm); radial distance from z-axis (mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_zbyr=new TH2D("HCal photon hit locations", "HCal photon hit locations;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_electron_zbyr=new TH2D("HCal electron hit locations", "HCal electron hit locations;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_neutron_zbyr=new TH2D("HCal neutron hit locations", "HCal neutron hit locations;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_other_zbyr=new TH2D("HCal other particle hit locations", "HCal other particle hit locations;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_unmatched_zbyr=new TH2D("HCal unmatched hit locations", "HCal unmatched hit locations;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_energy=new TH1F("HCal Photon hit Energies", "HCal Photon hit Energies;Energy(MeV);Count", 4000, 0, 4000);
        h_HCalhit_getTime=new TH1F("Creation time of particles causing HCal hits", "Creation time of particles causing HCal hits;Time(ns);Number of particles created", 500, 0, 500);
        h_HCalhit_getTime_nucleons=new TH1F("Creation time of nucleons causing HCal hits", "nucleons causing HCal hits;Time(ns);Number of Nucleons created", 500, 0, 500);
        h_HCalhit_nucleon_time_vs_energy=new TH2D("Nucleon time vs energy", "Nucleon time vs energy;Creation Time(ns);Energy(MeV)",100,0,500,250,0,4000); // 5ns time resolution, 16MeV resolution
        h_E_cal_summed_energy=new TH1F("E_cal_summed_energy","E_cal_summed_energy;Energy(MeV)(10MeV bin width);Count",800,0,8000);
        h_total_particles=new TH1F("total_particles","total_particles;Number of particles per event;Count",50,0,50);
        h_particle_energy=new TH1F("matched_particle_energy","matched_particle_energy;Energy(MeV)(5MeV bin width);Count",800,0,4000);
        h_hcal_hit_time_all=new TH1F("HCal_hit_time_all","HCal_hit_time_all;time(ns)(5ns bin width);Count",100,0,500);
        h_hit_time_creation_time_diff=new TH1F("hit_time_creation_time_diff","hit_time_creation_time_diff;time(ns)(2ns bin width);Count",100,0,200);
        h_part_hCalhit_tdif_less15_PE=new TH1F("part_hCalhit_tdif_less15_PE","part_hCalhit_tdif_less15_PE;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_less15_position=new TH2D("part_hCalhit_tdif_less15_position","part_hCalhit_tdif_less15_position;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_part_hCalhit_tdif_great40_PE=new TH1F("part_hCalhit_tdif_great40_PE","part_hCalhit_tdif_great40_PE;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_great40_position=new TH2D("part_hCalhit_tdif_great40_position","part_hCalhit_tdif_great40_position",80,0,3200,112,0,4500);
        h_hCalhit_time_less15_PE=new TH1F("hCalhit_time_less15_PE","hCalhit_time_less15_PE;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_less15_position=new TH2D("hCalhit_time_less15_position","hCalhit_time_less15_position",80,0,3200,112,0,4500);
        h_hCalhit_time_great40_PE=new TH1F("hCalhit_time_great40_PE","hCalhit_time_great40_PE;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_great40_position=new TH2D("hCalhit_time_great40_position","hCalhit_time_great40_position",80,0,3200,112,0,4500);
        h_hcal_hits_all_PEs=new TH1F("hcal_hits_all_PEs","hcal_hits_all_PEs;Photoelectrons(PEs);Count",200,0,200);
        h_hcal_hits_max_PE_of_event=new TH1F("h_hcal_hits_max_PE_of_event","h_hcal_hits_max_PE_of_event;Photoelectrons(PEs);Count",500,0,500);
    
        //SD region 1  (normal range, within ~1SD)
        h_PDGIDs_SD[0]=new TH1F("PDG_IDs_SD1","PDG_IDs_SD1",10000,-5000,5000);
        h_ZdepthofHCalHit_SD[0]=new TH1F("Z depth of HCal hits_SD1(10mm bins)", "Z depth of HCal hits_SD1(10mm bins)", 320, 0, 3200);
        h_ParticleHit_Distance_SD[0]=new TH1F("Distance between sim particle and HCal hits_SD1(5mm bins)", "Distance between sim particle and HCal hits_SD1(5mm bins)", 400, 0, 2000);
        h_HCalhit_zbyr_SD[0]=new TH2D("HCal hit locations_SD1", "HCal hit locations_SD1;Z depth (mm); radial distance from z-axis (mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_zbyr_SD[0]=new TH2D("HCal photon hit locations_SD1", "HCal photon hit locations_SD1;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_electron_zbyr_SD[0]=new TH2D("HCal electron hit locations_SD1", "HCal electron hit locations_SD1;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_neutron_zbyr_SD[0]=new TH2D("HCal neutron hit locations_SD1", "HCal neutron hit locations_SD1;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_other_zbyr_SD[0]=new TH2D("HCal other particle hit locations_SD1", "HCal other particle hit locations_SD1;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_unmatched_zbyr_SD[0]=new TH2D("HCal unmatched hit locations_SD1", "HCal unmatched hit locations_SD1;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_energy_SD[0]=new TH1F("HCal Photon hit Energies_SD1", "HCal Photon hit Energies_SD1;Energy(MeV);Count", 4000, 0, 4000);
        h_HCalhit_getTime_SD[0]=new TH1F("Creation time of particles causing HCal hits_SD1", "Creation time of particles causing HCal hits_SD1;Time(ns);Number of particles created", 500, 0, 500);
        h_HCalhit_getTime_nucleons_SD[0]=new TH1F("Creation time of nucleons causing HCal hits_SD1", "nucleons causing HCal hits_SD1;Time(ns);Number of Nucleons created", 500, 0, 500);
        h_HCalhit_nucleon_time_vs_energy_SD[0]=new TH2D("Nucleon time vs energy_SD1", "Nucleon time vs energy_SD1;Creation Time(ns);Energy(MeV)",100,0,500,250,0,4000); // 5ns time resolution, 16MeV resolution
        h_E_cal_summed_energy_SD[0]=new TH1F("E_cal_summed_energy_SD1","E_cal_summed_energy_SD1;Energy(MeV)(10MeV bin width);Count",800,0,8000);//10MeV bins
        h_total_particles_SD[0]=new TH1F("total_particles_SD1","total_particles_SD1;Number of particles per event;Count",50,0,50);
        h_particle_energy_SD[0]=new TH1F("matched_particle_energy_SD1","matched_particle_energy_SD1;Energy(MeV)(5MeV bin width);Count",800,0,4000);
        h_hcal_hit_time_all_SD[0]=new TH1F("HCal_hit_time_all_SD1","HCal_hit_time_all_SD1;time(ns)(5ns bin width);Count",100,0,500);
        h_hit_time_creation_time_diff_SD[0]=new TH1F("hit_time_creation_time_diff_SD1","hit_time_creation_time_diff_SD1;time(ns)(2ns bin width);Count",100,0,200);
        h_part_hCalhit_tdif_less15_PE_SD[0]=new TH1F("part_hCalhit_tdif_less15_PE_SD1","part_hCalhit_tdif_less15_PE_SD1;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_less15_position_SD[0]=new TH2D("part_hCalhit_tdif_less15_position_SD1","part_hCalhit_tdif_less15_position_SD1;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_part_hCalhit_tdif_great40_PE_SD[0]=new TH1F("part_hCalhit_tdif_great40_PE_SD1","part_hCalhit_tdif_great40_PE_SD1;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_great40_position_SD[0]=new TH2D("part_hCalhit_tdif_great40_position_SD1","part_hCalhit_tdif_great40_position_SD1",80,0,3200,112,0,4500);
        h_hCalhit_time_less15_PE_SD[0]=new TH1F("hCalhit_time_less15_PE_SD1","hCalhit_time_less15_PE_SD1;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_less15_position_SD[0]=new TH2D("hCalhit_time_less15_position_SD1","hCalhit_time_less15_position_SD1",80,0,3200,112,0,4500);
        h_hCalhit_time_great40_PE_SD[0]=new TH1F("hCalhit_time_great40_PE_SD1","hCalhit_time_great40_PE_SD1;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_great40_position_SD[0]=new TH2D("hCalhit_time_great40_position_SD1","hCalhit_time_great40_position_SD1",80,0,3200,112,0,4500);
        h_hcal_hits_all_PEs_SD[0]=new TH1F("hcal_hits_all_PEs_SD1","hcal_hits_all_PEs_SD1;Photoelectrons(PEs);Count",200,0,200);
        h_hcal_hits_max_PE_of_event_SD[0]=new TH1F("h_hcal_hits_max_PE_of_event_SD1","h_hcal_hits_max_PE_of_event_SD1;Photoelectrons(PEs);Count",500,0,500);
    
        //SD region 2 (high range, within than +1SD to +2SD)
        h_PDGIDs_SD[1]=new TH1F("PDG_IDs_SD2","PDG_IDs_SD2",10000,-5000,5000);
        h_ZdepthofHCalHit_SD[1]=new TH1F("Z depth of HCal hits_SD2(10mm bins)", "Z depth of HCal hits_SD2(10mm bins)", 320, 0, 3200);
        h_ParticleHit_Distance_SD[1]=new TH1F("Distance between sim particle and HCal hits_SD2(5mm bins)", "Distance between sim particle and HCal hits_SD2(5mm bins)", 400, 0, 2000);
        h_HCalhit_zbyr_SD[1]=new TH2D("HCal hit locations_SD2", "HCal hit locations_SD2;Z depth (mm); radial distance from z-axis (mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_zbyr_SD[1]=new TH2D("HCal photon hit locations_SD2", "HCal photon hit locations_SD2;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_electron_zbyr_SD[1]=new TH2D("HCal electron hit locations_SD2", "HCal electron hit locations_SD2;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_neutron_zbyr_SD[1]=new TH2D("HCal neutron hit locations_SD2", "HCal neutron hit locations_SD2;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_other_zbyr_SD[1]=new TH2D("HCal other particle hit locations_SD2", "HCal other particle hit locations_SD2;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_unmatched_zbyr_SD[1]=new TH2D("HCal unmatched hit locations_SD2", "HCal unmatched hit locations_SD2;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_energy_SD[1]=new TH1F("HCal Photon hit Energies_SD2", "HCal Photon hit Energies_SD2;Energy(MeV);Count", 4000, 0, 4000);
        h_HCalhit_getTime_SD[1]=new TH1F("Creation time of particles causing HCal hits_SD2", "Creation time of particles causing HCal hits_SD2;Time(ns);Number of particles created", 500, 0, 500);
        h_HCalhit_getTime_nucleons_SD[1]=new TH1F("Creation time of nucleons causing HCal hits_SD2", "nucleons causing HCal hits_SD2;Time(ns);Number of Nucleons created", 500, 0, 500);
        h_HCalhit_nucleon_time_vs_energy_SD[1]=new TH2D("Nucleon time vs energy_SD2", "Nucleon time vs energy_SD2;Creation Time(ns);Energy(MeV)",100,0,500,250,0,4000); // 5ns time resolution, 16MeV resolution
        h_E_cal_summed_energy_SD[1]=new TH1F("E_cal_summed_energy_SD2","E_cal_summed_energy_SD2;Energy(MeV)(10MeV bin width);Count",800,0,8000);//10MeV bins
        h_total_particles_SD[1]=new TH1F("total_particles_SD2","total_particles_SD2;Number of particles per event;Count",50,0,50);
        h_particle_energy_SD[1]=new TH1F("matched_particle_energy_SD2","matched_particle_energy_SD2;Energy(MeV)(5MeV bin width);Count",800,0,4000);
        h_hcal_hit_time_all_SD[1]=new TH1F("HCal_hit_time_all_SD2","HCal_hit_time_all_SD2;time(ns)(5ns bin width);Count",100,0,500);
        h_hit_time_creation_time_diff_SD[1]=new TH1F("hit_time_creation_time_diff_SD2","hit_time_creation_time_diff_SD2;time(ns)(2ns bin width);Count",100,0,200);
        h_part_hCalhit_tdif_less15_PE_SD[1]=new TH1F("part_hCalhit_tdif_less15_PE_SD2","part_hCalhit_tdif_less15_PE_SD2;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_less15_position_SD[1]=new TH2D("part_hCalhit_tdif_less15_position_SD2","part_hCalhit_tdif_less15_position_SD2;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_part_hCalhit_tdif_great40_PE_SD[1]=new TH1F("part_hCalhit_tdif_great40_PE_SD2","part_hCalhit_tdif_great40_PE_SD2;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_great40_position_SD[1]=new TH2D("part_hCalhit_tdif_great40_position_SD2","part_hCalhit_tdif_great40_position_SD2",80,0,3200,112,0,4500);
        h_hCalhit_time_less15_PE_SD[1]=new TH1F("hCalhit_time_less15_PE_SD2","hCalhit_time_less15_PE_SD2;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_less15_position_SD[1]=new TH2D("hCalhit_time_less15_position_SD2","hCalhit_time_less15_position_SD2",80,0,3200,112,0,4500);
        h_hCalhit_time_great40_PE_SD[1]=new TH1F("hCalhit_time_great40_PE_SD2","hCalhit_time_great40_PE_SD2;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_great40_position_SD[1]=new TH2D("hCalhit_time_great40_position_SD2","hCalhit_time_great40_position_SD2",80,0,3200,112,0,4500);
        h_hcal_hits_all_PEs_SD[1]=new TH1F("hcal_hits_all_PEs_SD2","hcal_hits_all_PEs_SD2;Photoelectrons(PEs);Count",200,0,200);
        h_hcal_hits_max_PE_of_event_SD[1]=new TH1F("h_hcal_hits_max_PE_of_event_SD2","h_hcal_hits_max_PE_of_event_SD2;Photoelectrons(PEs);Count",500,0,500);
    
        //SD region 3 (low range, within -1SD to -2SD)
        h_PDGIDs_SD[2]=new TH1F("PDG_IDs_SD3","PDG_IDs_SD3",10000,-5000,5000);
        h_ZdepthofHCalHit_SD[2]=new TH1F("Z depth of HCal hits_SD3(10mm bins)", "Z depth of HCal hits_SD3(10mm bins)", 320, 0, 3200);
        h_ParticleHit_Distance_SD[2]=new TH1F("Distance between sim particle and HCal hits_SD3(5mm bins)", "Distance between sim particle and HCal hits_SD3(5mm bins)", 400, 0, 2000);
        h_HCalhit_zbyr_SD[2]=new TH2D("HCal hit locations_SD3", "HCal hit locations_SD3;Z depth (mm); radial distance from z-axis (mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_zbyr_SD[2]=new TH2D("HCal photon hit locations_SD3", "HCal photon hit locations_SD3;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_electron_zbyr_SD[2]=new TH2D("HCal electron hit locations_SD3", "HCal electron hit locations_SD3;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_neutron_zbyr_SD[2]=new TH2D("HCal neutron hit locations_SD3", "HCal neutron hit locations_SD3;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_other_zbyr_SD[2]=new TH2D("HCal other particle hit locations_SD3", "HCal other particle hit locations_SD3;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_unmatched_zbyr_SD[2]=new TH2D("HCal unmatched hit locations_SD3", "HCal unmatched hit locations_SD3;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_energy_SD[2]=new TH1F("HCal Photon hit Energies_SD3", "HCal Photon hit Energies_SD3;Energy(MeV);Count", 4000, 0, 4000);
        h_HCalhit_getTime_SD[2]=new TH1F("Creation time of particles causing HCal hits_SD3", "Creation time of particles causing HCal hits_SD3;Time(ns);Number of particles created", 500, 0, 500);
        h_HCalhit_getTime_nucleons_SD[2]=new TH1F("Creation time of nucleons causing HCal hits_SD3", "nucleons causing HCal hits_SD3;Time(ns);Number of Nucleons created", 500, 0, 500);
        h_HCalhit_nucleon_time_vs_energy_SD[2]=new TH2D("Nucleon time vs energy_SD3", "Nucleon time vs energy_SD3;Creation Time(ns);Energy(MeV)",100,0,500,250,0,4000); // 5ns time resolution, 16MeV resolution
        h_E_cal_summed_energy_SD[2]=new TH1F("E_cal_summed_energy_SD3","E_cal_summed_energy_SD3;Energy(MeV)(10MeV bin width);Count",800,0,8000);//10MeV bins
        h_total_particles_SD[2]=new TH1F("total_particles_SD3","total_particles_SD3;Number of particles per event;Count",50,0,50);
        h_particle_energy_SD[2]=new TH1F("matched_particle_energy_SD3","matched_particle_energy_SD3;Energy(MeV)(5MeV bin width);Count",800,0,4000);
        h_hcal_hit_time_all_SD[2]=new TH1F("HCal_hit_time_all_SD3","HCal_hit_time_all_SD3;time(ns)(5ns bin width);Count",100,0,500);
        h_hit_time_creation_time_diff_SD[2]=new TH1F("hit_time_creation_time_diff_SD3","hit_time_creation_time_diff_SD3;time(ns)(2ns bin width);Count",100,0,200);
        h_part_hCalhit_tdif_less15_PE_SD[2]=new TH1F("part_hCalhit_tdif_less15_PE_SD3","part_hCalhit_tdif_less15_PE_SD3;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_less15_position_SD[2]=new TH2D("part_hCalhit_tdif_less15_position_SD3","part_hCalhit_tdif_less15_position_SD3;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_part_hCalhit_tdif_great40_PE_SD[2]=new TH1F("part_hCalhit_tdif_great40_PE_SD3","part_hCalhit_tdif_great40_PE_SD3;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_great40_position_SD[2]=new TH2D("part_hCalhit_tdif_great40_position_SD3","part_hCalhit_tdif_great40_position_SD3",80,0,3200,112,0,4500);
        h_hCalhit_time_less15_PE_SD[2]=new TH1F("hCalhit_time_less15_PE_SD3","hCalhit_time_less15_PE_SD3;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_less15_position_SD[2]=new TH2D("hCalhit_time_less15_position_SD3","hCalhit_time_less15_position_SD3",80,0,3200,112,0,4500);
        h_hCalhit_time_great40_PE_SD[2]=new TH1F("hCalhit_time_great40_PE_SD3","hCalhit_time_great40_PE_SD3;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_great40_position_SD[2]=new TH2D("hCalhit_time_great40_position_SD3","hCalhit_time_great40_position_SD3",80,0,3200,112,0,4500);
        h_hcal_hits_all_PEs_SD[2]=new TH1F("hcal_hits_all_PEs_SD3","hcal_hits_all_PEs_SD3;Photoelectrons(PEs);Count",200,0,200);
        h_hcal_hits_max_PE_of_event_SD[2]=new TH1F("h_hcal_hits_max_PE_of_event_SD3","h_hcal_hits_max_PE_of_event_SD3;Photoelectrons(PEs);Count",500,0,500);
    
        //SD region 4 (higher range, higher than 2SD)
        h_PDGIDs_SD[3]=new TH1F("PDG_IDs_SD4","PDG_IDs_SD4",10000,-5000,5000);
        h_ZdepthofHCalHit_SD[3]=new TH1F("Z depth of HCal hits_SD4(10mm bins)", "Z depth of HCal hits_SD4(10mm bins)", 320, 0, 3200);
        h_ParticleHit_Distance_SD[3]=new TH1F("Distance between sim particle and HCal hits_SD4(5mm bins)", "Distance between sim particle and HCal hits_SD4(5mm bins)", 400, 0, 2000);
        h_HCalhit_zbyr_SD[3]=new TH2D("HCal hit locations_SD4", "HCal hit locations_SD4;Z depth (mm); radial distance from z-axis (mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_zbyr_SD[3]=new TH2D("HCal photon hit locations_SD4", "HCal photon hit locations_SD4;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_electron_zbyr_SD[3]=new TH2D("HCal electron hit locations_SD4", "HCal electron hit locations_SD4;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_neutron_zbyr_SD[3]=new TH2D("HCal neutron hit locations_SD4", "HCal neutron hit locations_SD4;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_other_zbyr_SD[3]=new TH2D("HCal other particle hit locations_SD4", "HCal other particle hit locations_SD4;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_unmatched_zbyr_SD[3]=new TH2D("HCal unmatched hit locations_SD4", "HCal unmatched hit locations_SD4;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_energy_SD[3]=new TH1F("HCal Photon hit Energies_SD4", "HCal Photon hit Energies_SD4;Energy(MeV);Count", 4000, 0, 4000);
        h_HCalhit_getTime_SD[3]=new TH1F("Creation time of particles causing HCal hits_SD4", "Creation time of particles causing HCal hits_SD4;Time(ns);Number of particles created", 500, 0, 500);
        h_HCalhit_getTime_nucleons_SD[3]=new TH1F("Creation time of nucleons causing HCal hits_SD4", "nucleons causing HCal hits_SD4;Time(ns);Number of Nucleons created", 500, 0, 500);
        h_HCalhit_nucleon_time_vs_energy_SD[3]=new TH2D("Nucleon time vs energy_SD4", "Nucleon time vs energy_SD4;Creation Time(ns);Energy(MeV)",100,0,500,250,0,4000); // 5ns time resolution, 16MeV resolution
        h_E_cal_summed_energy_SD[3]=new TH1F("E_cal_summed_energy_SD4","E_cal_summed_energy_SD4;Energy(MeV)(10MeV bin width);Count",800,0,8000);//10MeV bins
        h_total_particles_SD[3]=new TH1F("total_particles_SD4","total_particles_SD4;Number of particles per event;Count",50,0,50);
        h_particle_energy_SD[3]=new TH1F("matched_particle_energy_SD4","matched_particle_energy_SD4;Energy(MeV)(5MeV bin width);Count",800,0,4000);
        h_hcal_hit_time_all_SD[3]=new TH1F("HCal_hit_time_all_SD4","HCal_hit_time_all_SD4;time(ns)(5ns bin width);Count",100,0,500);
        h_hit_time_creation_time_diff_SD[3]=new TH1F("hit_time_creation_time_diff_SD4","hit_time_creation_time_diff_SD4;time(ns)(2ns bin width);Count",100,0,200);
        h_part_hCalhit_tdif_less15_PE_SD[3]=new TH1F("part_hCalhit_tdif_less15_PE_SD4","part_hCalhit_tdif_less15_PE_SD4;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_less15_position_SD[3]=new TH2D("part_hCalhit_tdif_less15_position_SD4","part_hCalhit_tdif_less15_position_SD4;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_part_hCalhit_tdif_great40_PE_SD[3]=new TH1F("part_hCalhit_tdif_great40_PE_SD4","part_hCalhit_tdif_great40_PE_SD4;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_great40_position_SD[3]=new TH2D("part_hCalhit_tdif_great40_position_SD4","part_hCalhit_tdif_great40_position_SD4",80,0,3200,112,0,4500);
        h_hCalhit_time_less15_PE_SD[3]=new TH1F("hCalhit_time_less15_PE_SD4","hCalhit_time_less15_PE_SD4;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_less15_position_SD[3]=new TH2D("hCalhit_time_less15_position_SD4","hCalhit_time_less15_position_SD4",80,0,3200,112,0,4500);
        h_hCalhit_time_great40_PE_SD[3]=new TH1F("hCalhit_time_great40_PE_SD4","hCalhit_time_great40_PE_SD4;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_great40_position_SD[3]=new TH2D("hCalhit_time_great40_position_SD4","hCalhit_time_great40_position_SD4",80,0,3200,112,0,4500);
        h_hcal_hits_all_PEs_SD[3]=new TH1F("hcal_hits_all_PEs_SD4","hcal_hits_all_PEs_SD4;Photoelectrons(PEs);Count",200,0,200);
        h_hcal_hits_max_PE_of_event_SD[3]=new TH1F("h_hcal_hits_max_PE_of_event_SD4","h_hcal_hits_max_PE_of_event_SD4;Photoelectrons(PEs);Count",500,0,500);
    
        //SD region 5 (lower range, within -2SD to -3SD)
        h_PDGIDs_SD[4]=new TH1F("PDG_IDs_SD5","PDG_IDs_SD5",10000,-5000,5000);
        h_ZdepthofHCalHit_SD[4]=new TH1F("Z depth of HCal hits_SD5(10mm bins)", "Z depth of HCal hits_SD5(10mm bins)", 320, 0, 3200);
        h_ParticleHit_Distance_SD[4]=new TH1F("Distance between sim particle and HCal hits_SD5(5mm bins)", "Distance between sim particle and HCal hits_SD5(5mm bins)", 400, 0, 2000);
        h_HCalhit_zbyr_SD[4]=new TH2D("HCal hit locations_SD5", "HCal hit locations_SD5;Z depth (mm); radial distance from z-axis (mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_zbyr_SD[4]=new TH2D("HCal photon hit locations_SD5", "HCal photon hit locations_SD5;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_electron_zbyr_SD[4]=new TH2D("HCal electron hit locations_SD5", "HCal electron hit locations_SD5;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_neutron_zbyr_SD[4]=new TH2D("HCal neutron hit locations_SD5", "HCal neutron hit locations_SD5;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_other_zbyr_SD[4]=new TH2D("HCal other particle hit locations_SD5", "HCal other particle hit locations_SD5;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_unmatched_zbyr_SD[4]=new TH2D("HCal unmatched hit locations_SD5", "HCal unmatched hit locations_SD5;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_energy_SD[4]=new TH1F("HCal Photon hit Energies_SD5", "HCal Photon hit Energies_SD5;Energy(MeV);Count", 4000, 0, 4000);
        h_HCalhit_getTime_SD[4]=new TH1F("Creation time of particles causing HCal hits_SD5", "Creation time of particles causing HCal hits_SD5;Time(ns);Number of particles created", 500, 0, 500);
        h_HCalhit_getTime_nucleons_SD[4]=new TH1F("Creation time of nucleons causing HCal hits_SD5", "nucleons causing HCal hits_SD5;Time(ns);Number of Nucleons created", 500, 0, 500);
        h_HCalhit_nucleon_time_vs_energy_SD[4]=new TH2D("Nucleon time vs energy_SD5", "Nucleon time vs energy_SD5;Creation Time(ns);Energy(MeV)",100,0,500,250,0,4000); // 5ns time resolution, 16MeV resolution
        h_E_cal_summed_energy_SD[4]=new TH1F("E_cal_summed_energy_SD5","E_cal_summed_energy_SD5;Energy(MeV)(10MeV bin width);Count",800,0,8000);//10MeV bins
        h_total_particles_SD[4]=new TH1F("total_particles_SD5","total_particles_SD5;Number of particles per event;Count",50,0,50);
        h_particle_energy_SD[4]=new TH1F("matched_particle_energy_SD5","matched_particle_energy_SD5;Energy(MeV)(5MeV bin width);Count",800,0,4000);
        h_hcal_hit_time_all_SD[4]=new TH1F("HCal_hit_time_all_SD5","HCal_hit_time_all_SD5;time(ns)(5ns bin width);Count",100,0,500);
        h_hit_time_creation_time_diff_SD[4]=new TH1F("hit_time_creation_time_diff_SD5","hit_time_creation_time_diff_SD5;time(ns)(2ns bin width);Count",100,0,200);
        h_part_hCalhit_tdif_less15_PE_SD[4]=new TH1F("part_hCalhit_tdif_less15_PE_SD5","part_hCalhit_tdif_less15_PE_SD5;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_less15_position_SD[4]=new TH2D("part_hCalhit_tdif_less15_position_SD5","part_hCalhit_tdif_less15_position_SD5;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_part_hCalhit_tdif_great40_PE_SD[4]=new TH1F("part_hCalhit_tdif_great40_PE_SD5","part_hCalhit_tdif_great40_PE_SD5;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_great40_position_SD[4]=new TH2D("part_hCalhit_tdif_great40_position_SD5","part_hCalhit_tdif_great40_position_SD5",80,0,3200,112,0,4500);
        h_hCalhit_time_less15_PE_SD[4]=new TH1F("hCalhit_time_less15_PE_SD5","hCalhit_time_less15_PE_SD5;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_less15_position_SD[4]=new TH2D("hCalhit_time_less15_position_SD5","hCalhit_time_less15_position_SD5",80,0,3200,112,0,4500);
        h_hCalhit_time_great40_PE_SD[4]=new TH1F("hCalhit_time_great40_PE_SD5","hCalhit_time_great40_PE_SD5;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_great40_position_SD[4]=new TH2D("hCalhit_time_great40_position_SD5","hCalhit_time_great40_position_SD5",80,0,3200,112,0,4500);
        h_hcal_hits_all_PEs_SD[4]=new TH1F("hcal_hits_all_PEs_SD5","hcal_hits_all_PEs_SD5;Photoelectrons(PEs);Count",200,0,200);
        h_hcal_hits_max_PE_of_event_SD[4]=new TH1F("h_hcal_hits_max_PE_of_event_SD5","h_hcal_hits_max_PE_of_event_SD5;Photoelectrons(PEs);Count",500,0,500);
    
        //SD region 6 (very low range, within -3SD to -4SD)
        h_PDGIDs_SD[5]=new TH1F("PDG_IDs_SD6","PDG_IDs_SD6",10000,-5000,5000);
        h_ZdepthofHCalHit_SD[5]=new TH1F("Z depth of HCal hits_SD6(10mm bins)", "Z depth of HCal hits_SD6(10mm bins)", 320, 0, 3200);
        h_ParticleHit_Distance_SD[5]=new TH1F("Distance between sim particle and HCal hits_SD6(5mm bins)", "Distance between sim particle and HCal hits_SD6(5mm bins)", 400, 0, 2000);
        h_HCalhit_zbyr_SD[5]=new TH2D("HCal hit locations_SD6", "HCal hit locations_SD6;Z depth (mm); radial distance from z-axis (mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_zbyr_SD[5]=new TH2D("HCal photon hit locations_SD6", "HCal photon hit locations_SD6;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_electron_zbyr_SD[5]=new TH2D("HCal electron hit locations_SD6", "HCal electron hit locations_SD6;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_neutron_zbyr_SD[5]=new TH2D("HCal neutron hit locations_SD6", "HCal neutron hit locations_SD6;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_other_zbyr_SD[5]=new TH2D("HCal other particle hit locations_SD6", "HCal other particle hit locations_SD6;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_unmatched_zbyr_SD[5]=new TH2D("HCal unmatched hit locations_SD6", "HCal unmatched hit locations_SD6;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_energy_SD[5]=new TH1F("HCal Photon hit Energies_SD6", "HCal Photon hit Energies_SD6;Energy(MeV);Count", 4000, 0, 4000);
        h_HCalhit_getTime_SD[5]=new TH1F("Creation time of particles causing HCal hits_SD6", "Creation time of particles causing HCal hits_SD6;Time(ns);Number of particles created", 500, 0, 500);
        h_HCalhit_getTime_nucleons_SD[5]=new TH1F("Creation time of nucleons causing HCal hits_SD6", "nucleons causing HCal hits_SD6;Time(ns);Number of Nucleons created", 500, 0, 500);
        h_HCalhit_nucleon_time_vs_energy_SD[5]=new TH2D("Nucleon time vs energy_SD6", "Nucleon time vs energy_SD6;Creation Time(ns);Energy(MeV)",100,0,500,250,0,4000); // 5ns time resolution, 16MeV resolution
        h_E_cal_summed_energy_SD[5]=new TH1F("E_cal_summed_energy_SD6","E_cal_summed_energy_SD6;Energy(MeV)(10MeV bin width);Count",800,0,8000);//10MeV bins
        h_total_particles_SD[5]=new TH1F("total_particles_SD6","total_particles_SD6;Number of particles per event;Count",50,0,50);
        h_particle_energy_SD[5]=new TH1F("matched_particle_energy_SD6","matched_particle_energy_SD6;Energy(MeV)(5MeV bin width);Count",800,0,4000);
        h_hcal_hit_time_all_SD[5]=new TH1F("HCal_hit_time_all_SD6","HCal_hit_time_all_SD6;time(ns)(5ns bin width);Count",100,0,500);
        h_hit_time_creation_time_diff_SD[5]=new TH1F("hit_time_creation_time_diff_SD6","hit_time_creation_time_diff_SD6;time(ns)(2ns bin width);Count",100,0,200);
        h_part_hCalhit_tdif_less15_PE_SD[5]=new TH1F("part_hCalhit_tdif_less15_PE_SD6","part_hCalhit_tdif_less15_PE_SD6;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_less15_position_SD[5]=new TH2D("part_hCalhit_tdif_less15_position_SD6","part_hCalhit_tdif_less15_position_SD6;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_part_hCalhit_tdif_great40_PE_SD[5]=new TH1F("part_hCalhit_tdif_great40_PE_SD6","part_hCalhit_tdif_great40_PE_SD6;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_great40_position_SD[5]=new TH2D("part_hCalhit_tdif_great40_position_SD6","part_hCalhit_tdif_great40_position_SD6",80,0,3200,112,0,4500);
        h_hCalhit_time_less15_PE_SD[5]=new TH1F("hCalhit_time_less15_PE_SD6","hCalhit_time_less15_PE_SD6;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_less15_position_SD[5]=new TH2D("hCalhit_time_less15_position_SD6","hCalhit_time_less15_position_SD6",80,0,3200,112,0,4500);
        h_hCalhit_time_great40_PE_SD[5]=new TH1F("hCalhit_time_great40_PE_SD6","hCalhit_time_great40_PE_SD6;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_great40_position_SD[5]=new TH2D("hCalhit_time_great40_position_SD6","hCalhit_time_great40_position_SD6",80,0,3200,112,0,4500);
        h_hcal_hits_all_PEs_SD[5]=new TH1F("hcal_hits_all_PEs_SD6","hcal_hits_all_PEs_SD6;Photoelectrons(PEs);Count",200,0,200);
        h_hcal_hits_max_PE_of_event_SD[5]=new TH1F("h_hcal_hits_max_PE_of_event_SD6","h_hcal_hits_max_PE_of_event_SD6;Photoelectrons(PEs);Count",500,0,500);
    
        //SD region 7 (extremely low range, within -4SD to -5SD)
        h_PDGIDs_SD[6]=new TH1F("PDG_IDs_SD7","PDG_IDs_SD7",10000,-5000,5000);
        h_ZdepthofHCalHit_SD[6]=new TH1F("Z depth of HCal hits_SD7(10mm bins)", "Z depth of HCal hits_SD7(10mm bins)", 320, 0, 3200);
        h_ParticleHit_Distance_SD[6]=new TH1F("Distance between sim particle and HCal hits_SD7(5mm bins)", "Distance between sim particle and HCal hits_SD7(5mm bins)", 400, 0, 2000);
        h_HCalhit_zbyr_SD[6]=new TH2D("HCal hit locations_SD7", "HCal hit locations_SD7;Z depth (mm); radial distance from z-axis (mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_zbyr_SD[6]=new TH2D("HCal photon hit locations_SD7", "HCal photon hit locations_SD7;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_electron_zbyr_SD[6]=new TH2D("HCal electron hit locations_SD7", "HCal electron hit locations_SD7;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_neutron_zbyr_SD[6]=new TH2D("HCal neutron hit locations_SD7", "HCal neutron hit locations_SD7;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_other_zbyr_SD[6]=new TH2D("HCal other particle hit locations_SD7", "HCal other particle hit locations_SD7;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_unmatched_zbyr_SD[6]=new TH2D("HCal unmatched hit locations_SD7", "HCal unmatched hit locations_SD7;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_energy_SD[6]=new TH1F("HCal Photon hit Energies_SD7", "HCal Photon hit Energies_SD7;Energy(MeV);Count", 4000, 0, 4000);
        h_HCalhit_getTime_SD[6]=new TH1F("Creation time of particles causing HCal hits_SD7", "Creation time of particles causing HCal hits_SD7;Time(ns);Number of particles created", 500, 0, 500);
        h_HCalhit_getTime_nucleons_SD[6]=new TH1F("Creation time of nucleons causing HCal hits_SD7", "nucleons causing HCal hits_SD7;Time(ns);Number of Nucleons created", 500, 0, 500);
        h_HCalhit_nucleon_time_vs_energy_SD[6]=new TH2D("Nucleon time vs energy_SD7", "Nucleon time vs energy_SD7;Creation Time(ns);Energy(MeV)",100,0,500,250,0,4000); // 5ns time resolution, 16MeV resolution
        h_E_cal_summed_energy_SD[6]=new TH1F("E_cal_summed_energy_SD7","E_cal_summed_energy_SD7;Energy(MeV)(10MeV bin width);Count",800,0,8000);//10MeV bins
        h_total_particles_SD[6]=new TH1F("total_particles_SD7","total_particles_SD7;Number of particles per event;Count",50,0,50);
        h_particle_energy_SD[6]=new TH1F("matched_particle_energy_SD7","matched_particle_energy_SD7;Energy(MeV)(5MeV bin width);Count",800,0,4000);
        h_hcal_hit_time_all_SD[6]=new TH1F("HCal_hit_time_all_SD7","HCal_hit_time_all_SD7;time(ns)(5ns bin width);Count",100,0,500);
        h_hit_time_creation_time_diff_SD[6]=new TH1F("hit_time_creation_time_diff_SD7","hit_time_creation_time_diff_SD7;time(ns)(2ns bin width);Count",100,0,200);
        h_part_hCalhit_tdif_less15_PE_SD[6]=new TH1F("part_hCalhit_tdif_less15_PE_SD7","part_hCalhit_tdif_less15_PE_SD7;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_less15_position_SD[6]=new TH2D("part_hCalhit_tdif_less15_position_SD7","part_hCalhit_tdif_less15_position_SD7;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_part_hCalhit_tdif_great40_PE_SD[6]=new TH1F("part_hCalhit_tdif_great40_PE_SD7","part_hCalhit_tdif_great40_PE_SD7;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_great40_position_SD[6]=new TH2D("part_hCalhit_tdif_great40_position_SD7","part_hCalhit_tdif_great40_position_SD7",80,0,3200,112,0,4500);
        h_hCalhit_time_less15_PE_SD[6]=new TH1F("hCalhit_time_less15_PE_SD7","hCalhit_time_less15_PE_SD7;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_less15_position_SD[6]=new TH2D("hCalhit_time_less15_position_SD7","hCalhit_time_less15_position_SD7",80,0,3200,112,0,4500);
        h_hCalhit_time_great40_PE_SD[6]=new TH1F("hCalhit_time_great40_PE_SD7","hCalhit_time_great40_PE_SD7;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_great40_position_SD[6]=new TH2D("hCalhit_time_great40_position_SD7","hCalhit_time_great40_position_SD7",80,0,3200,112,0,4500);
        h_hcal_hits_all_PEs_SD[6]=new TH1F("hcal_hits_all_PEs_SD7","hcal_hits_all_PEs_SD7;Photoelectrons(PEs);Count",200,0,200);
        h_hcal_hits_max_PE_of_event_SD[6]=new TH1F("h_hcal_hits_max_PE_of_event_SD7","h_hcal_hits_max_PE_of_event_SD7;Photoelectrons(PEs);Count",500,0,500);
    
        //SD region 8 (super-duper low range, within -5SD to -6SD)
        h_PDGIDs_SD[7]=new TH1F("PDG_IDs_SD8","PDG_IDs_SD8",10000,-5000,5000);
        h_ZdepthofHCalHit_SD[7]=new TH1F("Z depth of HCal hits_SD8(10mm bins)", "Z depth of HCal hits_SD8(10mm bins)", 320, 0, 3200);
        h_ParticleHit_Distance_SD[7]=new TH1F("Distance between sim particle and HCal hits_SD8(5mm bins)", "Distance between sim particle and HCal hits_SD8(5mm bins)", 400, 0, 2000);
        h_HCalhit_zbyr_SD[7]=new TH2D("HCal hit locations_SD8", "HCal hit locations_SD8;Z depth (mm); radial distance from z-axis (mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_zbyr_SD[7]=new TH2D("HCal photon hit locations_SD8", "HCal photon hit locations_SD8;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_electron_zbyr_SD[7]=new TH2D("HCal electron hit locations_SD8", "HCal electron hit locations_SD8;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_neutron_zbyr_SD[7]=new TH2D("HCal neutron hit locations_SD8", "HCal neutron hit locations_SD8;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_other_zbyr_SD[7]=new TH2D("HCal other particle hit locations_SD8", "HCal other particle hit locations_SD8;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_unmatched_zbyr_SD[7]=new TH2D("HCal unmatched hit locations_SD8", "HCal unmatched hit locations_SD8;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_energy_SD[7]=new TH1F("HCal Photon hit Energies_SD8", "HCal Photon hit Energies_SD8;Energy(MeV);Count", 4000, 0, 4000);
        h_HCalhit_getTime_SD[7]=new TH1F("Creation time of particles causing HCal hits_SD8", "Creation time of particles causing HCal hits_SD8;Time(ns);Number of particles created", 500, 0, 500);
        h_HCalhit_getTime_nucleons_SD[7]=new TH1F("Creation time of nucleons causing HCal hits_SD8", "nucleons causing HCal hits_SD8;Time(ns);Number of Nucleons created", 500, 0, 500);
        h_HCalhit_nucleon_time_vs_energy_SD[7]=new TH2D("Nucleon time vs energy_SD8", "Nucleon time vs energy_SD8;Creation Time(ns);Energy(MeV)",100,0,500,250,0,4000); // 5ns time resolution, 16MeV resolution
        h_E_cal_summed_energy_SD[7]=new TH1F("E_cal_summed_energy_SD8","E_cal_summed_energy_SD8;Energy(MeV)(10MeV bin width);Count",800,0,8000);//10MeV bins
        h_total_particles_SD[7]=new TH1F("total_particles_SD8","total_particles_SD8;Number of particles per event;Count",50,0,50);
        h_particle_energy_SD[7]=new TH1F("matched_particle_energy_SD8","matched_particle_energy_SD8;Energy(MeV)(5MeV bin width);Count",800,0,4000);
        h_hcal_hit_time_all_SD[7]=new TH1F("HCal_hit_time_all_SD8","HCal_hit_time_all_SD8;time(ns)(5ns bin width);Count",100,0,500);
        h_hit_time_creation_time_diff_SD[7]=new TH1F("hit_time_creation_time_diff_SD8","hit_time_creation_time_diff_SD8;time(ns)(2ns bin width);Count",100,0,200);
        h_part_hCalhit_tdif_less15_PE_SD[7]=new TH1F("part_hCalhit_tdif_less15_PE_SD8","part_hCalhit_tdif_less15_PE_SD8;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_less15_position_SD[7]=new TH2D("part_hCalhit_tdif_less15_position_SD8","part_hCalhit_tdif_less15_position_SD8;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_part_hCalhit_tdif_great40_PE_SD[7]=new TH1F("part_hCalhit_tdif_great40_PE_SD8","part_hCalhit_tdif_great40_PE_SD8;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_great40_position_SD[7]=new TH2D("part_hCalhit_tdif_great40_position_SD8","part_hCalhit_tdif_great40_position_SD8",80,0,3200,112,0,4500);
        h_hCalhit_time_less15_PE_SD[7]=new TH1F("hCalhit_time_less15_PE_SD8","hCalhit_time_less15_PE_SD8;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_less15_position_SD[7]=new TH2D("hCalhit_time_less15_position_SD8","hCalhit_time_less15_position_SD8",80,0,3200,112,0,4500);
        h_hCalhit_time_great40_PE_SD[7]=new TH1F("hCalhit_time_great40_PE_SD8","hCalhit_time_great40_PE_SD8;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_great40_position_SD[7]=new TH2D("hCalhit_time_great40_position_SD8","hCalhit_time_great40_position_SD8",80,0,3200,112,0,4500);
        h_hcal_hits_all_PEs_SD[7]=new TH1F("hcal_hits_all_PEs_SD8","hcal_hits_all_PEs_SD8;Photoelectrons(PEs);Count",200,0,200);
        h_hcal_hits_max_PE_of_event_SD[7]=new TH1F("h_hcal_hits_max_PE_of_event_SD8","h_hcal_hits_max_PE_of_event_SD8;Photoelectrons(PEs);Count",500,0,500);
    
        //SD region 9 (mega-ultra-super-duper low range, within -6SD to -7SD)
        h_PDGIDs_SD[8]=new TH1F("PDG_IDs_SD9","PDG_IDs_SD9",10000,-5000,5000);
        h_ZdepthofHCalHit_SD[8]=new TH1F("Z depth of HCal hits_SD9(10mm bins)", "Z depth of HCal hits_SD9(10mm bins)", 320, 0, 3200);
        h_ParticleHit_Distance_SD[8]=new TH1F("Distance between sim particle and HCal hits_SD9(5mm bins)", "Distance between sim particle and HCal hits_SD9(5mm bins)", 400, 0, 2000);
        h_HCalhit_zbyr_SD[8]=new TH2D("HCal hit locations_SD9", "HCal hit locations_SD9;Z depth (mm); radial distance from z-axis (mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_zbyr_SD[8]=new TH2D("HCal photon hit locations_SD9", "HCal photon hit locations_SD9;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_electron_zbyr_SD[8]=new TH2D("HCal electron hit locations_SD9", "HCal electron hit locations_SD9;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_neutron_zbyr_SD[8]=new TH2D("HCal neutron hit locations_SD9", "HCal neutron hit locations_SD9;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_other_zbyr_SD[8]=new TH2D("HCal other particle hit locations_SD9", "HCal other particle hit locations_SD9;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_unmatched_zbyr_SD[8]=new TH2D("HCal unmatched hit locations_SD9", "HCal unmatched hit locations_SD9;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_HCalhit_photon_energy_SD[8]=new TH1F("HCal Photon hit Energies_SD9", "HCal Photon hit Energies_SD9;Energy(MeV);Count", 4000, 0, 4000);
        h_HCalhit_getTime_SD[8]=new TH1F("Creation time of particles causing HCal hits_SD9", "Creation time of particles causing HCal hits_SD9;Time(ns);Number of particles created", 500, 0, 500);
        h_HCalhit_getTime_nucleons_SD[8]=new TH1F("Creation time of nucleons causing HCal hits_SD9", "nucleons causing HCal hits_SD9;Time(ns);Number of Nucleons created", 500, 0, 500);
        h_HCalhit_nucleon_time_vs_energy_SD[8]=new TH2D("Nucleon time vs energy_SD9", "Nucleon time vs energy_SD9;Creation Time(ns);Energy(MeV)",100,0,500,250,0,4000); // 5ns time resolution, 16MeV resolution
        h_E_cal_summed_energy_SD[8]=new TH1F("E_cal_summed_energy_SD9","E_cal_summed_energy_SD9;Energy(MeV)(10MeV bin width);Count",800,0,8000);//10MeV bins
        h_total_particles_SD[8]=new TH1F("total_particles_SD9","total_particles_SD9;Number of particles per event;Count",50,0,50);
        h_particle_energy_SD[8]=new TH1F("matched_particle_energy_SD9","matched_particle_energy_SD9;Energy(MeV)(5MeV bin width);Count",800,0,4000);
        h_hcal_hit_time_all_SD[8]=new TH1F("HCal_hit_time_all_SD9","HCal_hit_time_all_SD9;time(ns)(5ns bin width);Count",100,0,500);
        h_hit_time_creation_time_diff_SD[8]=new TH1F("hit_time_creation_time_diff_SD9","hit_time_creation_time_diff_SD9;time(ns)(2ns bin width);Count",100,0,200);
        h_part_hCalhit_tdif_less15_PE_SD[8]=new TH1F("part_hCalhit_tdif_less15_PE_SD9","part_hCalhit_tdif_less15_PE_SD9;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_less15_position_SD[8]=new TH2D("part_hCalhit_tdif_less15_position_SD9","part_hCalhit_tdif_less15_position_SD9;Z depth(mm);radial distance from z-axis(mm)",80,0,3200,112,0,4500);
        h_part_hCalhit_tdif_great40_PE_SD[8]=new TH1F("part_hCalhit_tdif_great40_PE_SD9","part_hCalhit_tdif_great40_PE_SD9;Photoelectrons(PEs);Count",200,0,200);
        h_part_hCalhit_tdif_great40_position_SD[8]=new TH2D("part_hCalhit_tdif_great40_position_SD9","part_hCalhit_tdif_great40_position_SD9",80,0,3200,112,0,4500);
        h_hCalhit_time_less15_PE_SD[8]=new TH1F("hCalhit_time_less15_PE_SD9","hCalhit_time_less15_PE_SD9;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_less15_position_SD[8]=new TH2D("hCalhit_time_less15_position_SD9","hCalhit_time_less15_position_SD9",80,0,3200,112,0,4500);
        h_hCalhit_time_great40_PE_SD[8]=new TH1F("hCalhit_time_great40_PE_SD9","hCalhit_time_great40_PE_SD9;Photoelectrons(PEs);Count",200,0,200);
        h_hCalhit_time_great40_position_SD[8]=new TH2D("hCalhit_time_great40_position_SD9","hCalhit_time_great40_position_SD9",80,0,3200,112,0,4500);
        h_hcal_hits_all_PEs_SD[8]=new TH1F("hcal_hits_all_PEs_SD9","hcal_hits_all_PEs_SD9;Photoelectrons(PEs);Count",200,0,200);
        h_hcal_hits_max_PE_of_event_SD[8]=new TH1F("h_hcal_hits_max_PE_of_event_SD9","h_hcal_hits_max_PE_of_event_SD9;Photoelectrons(PEs);Count",500,0,500);

        return;
    } //onProcessStart
    
} //ldmx namespace

DECLARE_ANALYZER_NS(ldmx, HcalHitMatcher);
