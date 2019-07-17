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
            if (sP == lastP) continue; //skip repeats
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
        h_Ecal_SummedEnergy_SD[9]->Fill(e_cal_sum_energy);
        h_Ecal_SummedEnergy_SD[ecal_sumESD]->Fill(e_cal_sum_energy);
    
        h_NumParticles_SD[9]->Fill(filteredSimVec.size());
        h_NumParticles_SD[ecal_sumESD]->Fill(filteredSimVec.size());

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
                    
                    h_Particle_HitDistance_All_SD[9]->Fill(new_dist);
                    h_Particle_HitDistance_All_SD[ecal_sumESD]->Fill(new_dist);
        
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

                h_HcalHit_ZbyR_All_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                h_HcalHit_Z_SD[9]->Fill(hcalhit->getZ());
                h_HcalHit_Time_All_SD[9]->Fill(hcalhit->getTime());
                h_HcalHit_PE_All_SD[9]->Fill(hcalhit->getPE());
                
                h_HcalHit_ZbyR_All_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                h_HcalHit_Z_SD[ecal_sumESD]->Fill(hcalhit->getZ());
                h_HcalHit_Time_All_SD[ecal_sumESD]->Fill(hcalhit->getTime());
                h_HcalHit_PE_All_SD[ecal_sumESD]->Fill(hcalhit->getPE());
                
                if(hcalhit->getTime() < 15.0)  {
                    h_HcalHit_PE_TimeLess15_SD[9]->Fill(hcalhit->getPE());
                    h_HcalHit_PE_TimeLess15_SD[ecal_sumESD]->Fill(hcalhit->getPE());
                    h_HcalHit_ZbyR_TimeLess15_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                    h_HcalHit_ZbyR_TimeLess15_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                } else if(hcalhit->getTime() > 40.0)  {
                    h_HcalHit_PE_TimeGreat40_SD[9]->Fill(hcalhit->getPE());
                    h_HcalHit_PE_TimeGreat40_SD[ecal_sumESD]->Fill(hcalhit->getPE());
                    h_HcalHit_ZbyR_TimeGreat40_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                    h_HcalHit_ZbyR_TimeGreat40_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                }
                
                if(hcalhit->getPE() > max_PE_of_event)
                    max_PE_of_event=hcalhit->getPE();
                
                if( dist <= 150.0 ) {//must be 150mm or closer to confidentally match an HCal hit to a sim particle
                
                    numMatchedHits_++;

                    h_Particle_HitDistance_Matched_SD[9]->Fill(new_dist);
                    h_Particle_HitDistance_Matched_SD[ecal_sumESD]->Fill(new_dist);
        
                    h_HcalHit_Time_Matched_All_SD[9]->Fill(hcalhit->getTime());
                    h_HcalHit_Time_Matched_All_SD[ecal_sumESD]->Fill(hcalhit->getTime());
                    
                    h_Particle_PDGID_Matched_SD[9]->Fill(pdgID);
                    h_Particle_PDGID_Matched_SD[ecal_sumESD]->Fill(pdgID);
        
                    double part_hcalhit_timeDiff = (hcalhit->getTime()) - (filteredSimVec[simPartNum]->getSimParticle()->getTime());
                    
                    h_HcalHit_Time_Matched_Tdiff_SD[9]->Fill(part_hcalhit_timeDiff);
                    h_HcalHit_Time_Matched_Tdiff_SD[ecal_sumESD]->Fill(part_hcalhit_timeDiff);

                    if(part_hcalhit_timeDiff < 15.0)  {
                        h_HcalHit_PE_Matched_TdifLess15_SD[9]->Fill(hcalhit->getPE());
                        h_HcalHit_PE_Matched_TdifLess15_SD[ecal_sumESD]->Fill(hcalhit->getPE());
                        h_HcalHit_ZbyR_Matched_TdifLess15_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                        h_HcalHit_ZbyR_Matched_TdifLess15_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                    } else if(part_hcalhit_timeDiff > 40.0)  {
                        h_HcalHit_PE_Matched_TdifGreat40_SD[9]->Fill(hcalhit->getPE());
                        h_HcalHit_PE_Matched_TdifGreat40_SD[ecal_sumESD]->Fill(hcalhit->getPE());
                        h_HcalHit_ZbyR_Matched_TdifGreat40_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                        h_HcalHit_ZbyR_Matched_TdifGreat40_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                    }
                    
                    //protons or neutrons (nucleons) 
                    if( pdgID==2112 or pdgID==2212 ) { 
                        h_HcalHit_Time_Matched_Nucleons_SD[9]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getTime());
//                        h_HCalhit_nucleon_time_vs_energy_SD[9]->Fill(
//                                filteredSimVec[simPartNum]->getSimParticle()->getTime(),
//                                filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                        h_HcalHit_Time_Matched_Nucleons_SD[ecal_sumESD]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getTime());
//                        h_HCalhit_nucleon_time_vs_energy_SD[ecal_sumESD]->Fill(
//                                filteredSimVec[simPartNum]->getSimParticle()->getTime(),
//                                filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                    }
                    
                    h_Particle_Energy_Matched_SD[9]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                    h_Particle_Energy_Matched_SD[ecal_sumESD]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
        
                    switch(pdgID) {
                        case 11:
                            h_HcalHit_ZbyR_Matched_Electron_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            h_HcalHit_ZbyR_Matched_Electron_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            break;
                        case 22:
                            h_HcalHit_ZbyR_Matched_Photon_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            h_HcalHit_ZbyR_Matched_Photon_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
//                            h_HCalhit_photon_energy_SD[9]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
//                            h_HCalhit_photon_energy_SD[ecal_sumESD]->Fill(filteredSimVec[simPartNum]->getSimParticle()->getEnergy());
                            break;
                        case 2112:
                            h_HcalHit_ZbyR_Matched_Neutron_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                            h_HcalHit_ZbyR_Matched_Neutron_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            break;
                        default:
                            h_HcalHit_ZbyR_Matched_Other_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            h_HcalHit_ZbyR_Matched_Other_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            break;
                    }

                } else {

                    h_HcalHit_ZbyR_Unmatched_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist);
                    h_HcalHit_ZbyR_Unmatched_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);

                } //matched or unmatched
    
            } // if not a noise hit

        }//End loop over hcalhits array

        // maximum PE in hcal hits for the event
        h_EventMaxPE_SD[9]->Fill(max_PE_of_event);
        h_EventMaxPE_SD[ecal_sumESD]->Fill(max_PE_of_event);

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

            h_Ecal_SummedEnergy_SD[i]=new TH1D(
                    ("Ecal_SummedEnergy_SD_"+range_name).c_str(),
                    ("Ecal Summed Energy SD "+range_name+";Energy(MeV) (10MeV bin width);Count").c_str(),
                    800,0,8000);//10MeV bins

            h_NumParticles_SD[i]=new TH1D(
                    ("NumParticles_SD_"+range_name).c_str(),
                    ("Num Particles in Event SD "+range_name+";Number of particles per event;Count").c_str(),
                    50,0,50);

            h_EventMaxPE_SD[i]=new TH1D(
                    ("EventMaxPE_SD_"+range_name).c_str(),
                    ("Maximum PE for all Hits in Event SD "+range_name+";Photoelectrons(PEs);Count").c_str(),
                    500,0,500);

            h_Particle_PDGID_All_SD[i]=new TH1D(
                    ("Particle_PDGID_All_SD_"+range_name).c_str(),
                    ("PDG IDs SD "+range_name).c_str(),
                    10000,-5000,5000);

            h_Particle_PDGID_Matched_SD[i]=new TH1D(
                    ("Particle_PDGID_Matched_SD_"+range_name).c_str(),
                    ("PDG IDs when matched SD "+range_name).c_str(),
                    10000,-5000,5000);

            h_Particle_HitDistance_All_SD[i]=new TH1D(
                    ("Particle_HitDistance_All_SD_"+range_name).c_str(), 
                    ("Distance between SimParticle and HcalHit SD "+range_name+" (5mm bins)").c_str(), 
                    400, 0, 2000);
            
            h_Particle_HitDistance_Matched_SD[i]=new TH1D(
                    ("Particle_HitDistance_Matched_SD_"+range_name).c_str(), 
                    ("Distance between SimParticle and HcalHit when matched SD "+range_name+" (5mm bins)").c_str(), 
                    400, 0, 2000);
            
            h_Particle_Energy_All_SD[i]=new TH1D(
                    ("Particle_Energy_All_SD_"+range_name).c_str(),
                    ("All Particle Energies SD "+range_name+";Energy(MeV)(5MeV bin width);Count").c_str(),
                    800,0,4000);

            h_Particle_Energy_Matched_SD[i]=new TH1D(
                    ("Particle_Energy_Matched_SD_"+range_name).c_str(),
                    ("Matched Particle Energies SD "+range_name+";Energy(MeV)(5MeV bin width);Count").c_str(),
                    800,0,4000);

            h_HcalHit_Z_SD[i]=new TH1D(
                    ("HcalHit_Z_SD_"+range_name).c_str(),
                    ("Z depth of HCal hits SD "+range_name+" (10mm bins)").c_str(),
                    320, 0, 3200);

            h_HcalHit_ZbyR_All_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_All_SD_"+range_name).c_str(), 
                    ("All HcalHit locations SD "+range_name+";Z depth (mm); radial distance from z-axis (mm)").c_str(),
                    80,0,3200,112,0,4500);

            h_HcalHit_ZbyR_Unmatched_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Unmatched_SD_"+range_name).c_str(), 
                    ("Hcal unmatched hit locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            h_HcalHit_ZbyR_TimeLess15_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_TimeLess15_SD_"+range_name).c_str(),
                    ("HcalHits with Time < 15ns locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            h_HcalHit_ZbyR_TimeGreat40_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_TimeGreat40_SD_"+range_name).c_str(),
                    ("HcalHits with Time > 40ns locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            h_HcalHit_ZbyR_Matched_Photon_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Matched_Photon_SD_"+range_name).c_str(), 
                    ("Hcal photon hit locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            h_HcalHit_ZbyR_Matched_Electron_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Matched_Electron_SD_"+range_name).c_str(), 
                    ("Hcal electron hit locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            h_HcalHit_ZbyR_Matched_Neutron_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Matched_Neutron_SD_"+range_name).c_str(), 
                    ("Hcal neutron hit locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            h_HcalHit_ZbyR_Matched_Other_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Matched_Other_SD_"+range_name).c_str(), 
                    ("Hcal other particle hit locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            h_HcalHit_ZbyR_Matched_TdifLess15_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Matched_TdifLess15_SD_"+range_name).c_str(),
                    ("Matched HcalHit location with time diff < 15ns SD "+
                     range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            h_HcalHit_ZbyR_Matched_TdifGreat40_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Matched_TdifGreat40_SD_"+range_name).c_str(),
                    ("Matched HcalHit location with time diff > 40ns SD "+
                     range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            h_HcalHit_PE_All_SD[i]=new TH1D(
                    ("HcalHit_PE_All_SD_"+range_name).c_str(),
                    ("PEs of all HcalHits SD "+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            h_HcalHit_PE_TimeLess15_SD[i]=new TH1D(
                    ("HcalHit_PE_TimeLess15_SD_"+range_name).c_str(),
                    ("HcalHits with Time < 15ns SD "+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            h_HcalHit_PE_TimeGreat40_SD[i]=new TH1D(
                    ("HcalHit_PE_TimeGreat40_SD_"+range_name).c_str(),
                    ("HcalHits with Time > 40ns SD "+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            h_HcalHit_PE_Matched_TdifLess15_SD[i]=new TH1D(
                    ("HcalHit_PE_Matched_TdifLess15_SD_"+range_name).c_str(),
                    ("Matched HcalHit with Time diff < 15ns SD "+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            h_HcalHit_PE_Matched_TdifGreat40_SD[i]=new TH1D(
                    ("HcalHit_PE_Matched_TdifGreat40_SD_"+range_name).c_str(),
                    ("Matched HcalHit with Time diff > 40ns SD "+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            h_HcalHit_Time_All_SD[i]=new TH1D(
                    ("HcalHit_Time_All_SD_"+range_name).c_str(),
                    ("Time of All HcalHits SD "+range_name+";time(ns)(5ns bin width);Count").c_str(),
                    100,0,500);

            h_HcalHit_Time_Matched_All_SD[i]=new TH1D(
                    ("HcalHit_Time_Matched_All_SD_"+range_name).c_str(), 
                    ("Time of Matched HcalHits SD "+range_name+";Time(ns);Number of particles created").c_str(), 
                    500, 0, 500);

            h_HcalHit_Time_Matched_Nucleons_SD[i]=new TH1D(
                    ("HcalHit_Time_Matched_Nucleons_SD_"+range_name).c_str(),
                    ("Time of HcalHits Matched to Nucleons SD "+range_name+";Time(ns);Number of Nucleons created").c_str(), 
                    500, 0, 500);
            
            h_HcalHit_Time_Matched_Tdiff_SD[i]=new TH1D(
                    ("HcalHit_Time_Matched_Tdiff_SD_"+range_name).c_str(),
                    ("Time Diff between SimParticle and matched HcalHit SD "+range_name+";time(ns)(2ns bin width);Count").c_str(),
                    100,0,200);

//            h_HCalhit_photon_energy_SD[i]=new TH1D(
//                    ("HCal Photon hit Energies_SD_"+range_name).c_str(),
//                    ("HCal Photon hit Energies_SD_"+range_name+";Energy(MeV);Count").c_str(), 
//                    4000, 0, 4000);
//
//            h_HCalhit_nucleon_time_vs_energy_SD[i]=new TH2D(
//                    ("Nucleon time vs energy_SD_"+range_name).c_str(), 
//                    ("Nucleon time vs energy_SD_"+range_name+";Creation Time(ns);Energy(MeV)").c_str(),
//                    100,0,500,250,0,4000); // 5ns time resolution, 16MeV resolution

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
