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
        
        EcalHitColl_ = ps.getString( "EcalHitCollectionName" , "ecalDigis" );
        HcalHitColl_ = ps.getString( "HcalHitCollectionName" , "hcalDigis" );
        EcalScoringPlane_ = ps.getString( "EcalScoringPlaneHitsName" , "EcalScoringPlaneHits" ); 
        HcalScoringPlane_ = ps.getString( "HcalScoringPlaneHitsName" , "HcalScoringPlaneHits" );

        maxMatchDist_ = ps.getDouble( "MaximumMatchDistance" , 150.0 );

        return;
    }

    void HcalHitMatcher::analyze(const ldmx::Event& event) {

        numEvents_++;

        //----------This section obtains a list of sim particles that cross the ecal scoring plane---------->
        const TClonesArray* ecalScoringPlaneHits = event.getCollection( EcalScoringPlane_ );
        const TClonesArray* simParticles = event.getCollection("SimParticles"); // side-effect is to make TRefs all valid

        std::vector<ldmx::SimTrackerHit*> simVec;
    
        for (int i = 0; i < ecalScoringPlaneHits->GetEntriesFast(); i++ ) {
            ldmx::SimTrackerHit* ecalSPH = (ldmx::SimTrackerHit*)(ecalScoringPlaneHits->At(i));
            simVec.push_back(ecalSPH);
        }

        std::sort(simVec.begin(), simVec.end(), compSims);

        //SimParticles that cross Ecal Scoring Planes
        std::vector< ldmx::SimParticle* > simParticleCrossEcalSP;

        ldmx::SimParticle* lastP = 0; //sometimes multiple SP hits from same particle
        for ( std::vector<ldmx::SimTrackerHit*>::iterator it_simVec = simVec.begin();
              it_simVec != simVec.end(); ++it_simVec ) {

            ldmx::SimParticle* sP = (*it_simVec)->getSimParticle();
            if ( sP != lastP ) {
                //make sure there aren't any repeats
                lastP = sP;
                simParticleCrossEcalSP.push_back( sP );
            }
        }

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
    
        h_NumParticles_SD[9]->Fill(simParticleCrossEcalSP.size());
        h_NumParticles_SD[ecal_sumESD]->Fill(simParticleCrossEcalSP.size());

        //Go through all SimParticles that crossed ECAL Scoring Plane
        for ( const ldmx::SimParticle* simPart : simParticleCrossEcalSP ) {
            int pdgID = simPart->getPdgID();
            double energy = simPart->getEnergy();

            h_Particle_PDGID_All_SD[9]->Fill( pdgID );
            h_Particle_PDGID_All_SD[ecal_sumESD]->Fill( pdgID );

            h_Particle_Energy_All_SD[9]->Fill( energy );
            h_Particle_Energy_All_SD[ecal_sumESD]->Fill( energy );
        }

        //----This section matches HCal hits to sim particles and records results----->
        const TClonesArray* hcalHitColl = event.getCollection( HcalHitColl_ );

        float max_PE_of_event=0;
        for(int i=0; i < hcalHitColl->GetEntriesFast(); i++) { //Begin loop over hcalhits array
            ldmx::HcalHit* hcalhit = (ldmx::HcalHit*)(hcalHitColl->At(i));
            
            if ( ! hcalhit->getNoise() ) { //Only analyze non-noise hits

                numNonNoiseHits_++;

                //---- Bin HcalHit information that does not depend on mathcin -------------------->
                
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

                //---- Attempt to match this HcalHit to a SimParticle that cross the Ecal SP ----->

                //Iterate over all SimParticles that cross Ecal Scoring Plane to try to find a match
                //  for this HcalHit
                double new_dist=9999, dist=9998; //initializing distance variables
                const ldmx::SimParticle* matchedParticle = nullptr;
                for ( const ldmx::SimParticle* simPart : simParticleCrossEcalSP ) {

                    std::vector<double> simStart = simPart->getVertex();
                    std::vector<double> simEnd = simPart->getEndPoint();
        
                    TVector3 simStartT = TVector3(simStart[0], simStart[1], simStart[2]);
                    TVector3 simEndT = TVector3(simEnd[0], simEnd[1], simEnd[2]);
                    TVector3 hCalPoint = TVector3(hcalhit->getX(), hcalhit->getY(), hcalhit->getZ());
                    
                    new_dist = point_line_distance(simStartT, simEndT, hCalPoint);
                    
                    h_Particle_HitDistance_All_SD[9]->Fill(new_dist);
                    h_Particle_HitDistance_All_SD[ecal_sumESD]->Fill(new_dist);
        
                    if(simStart[2]<10.0 and simPart->getEnergy()>3000.0) {
                        //discarding original electron
                        //DO NOTHING (skip this sim particle)
                    } else if(new_dist < dist) {
                        dist = new_dist; //Distance to matched particle
                        matchedParticle = simPart;
                    }
                    
                } //iterate over sim particles to match one to current hcal hit

                //---- Bin HcalHit/SimParticle information for successfully matched hits --------->

                if( matchedParticle and dist <= maxMatchDist_ ) {
                
                    numMatchedHits_++;

                    h_Particle_HitDistance_Matched_SD[9]->Fill( dist );
                    h_Particle_HitDistance_Matched_SD[ecal_sumESD]->Fill( dist );
        
                    h_HcalHit_Time_Matched_All_SD[9]->Fill(hcalhit->getTime());
                    h_HcalHit_Time_Matched_All_SD[ecal_sumESD]->Fill(hcalhit->getTime());
                    
                    int pdgID = matchedParticle->getPdgID();

                    h_Particle_PDGID_Matched_SD[9]->Fill(pdgID);
                    h_Particle_PDGID_Matched_SD[ecal_sumESD]->Fill(pdgID);
        
                    double part_hcalhit_timeDiff = (hcalhit->getTime()) - (matchedParticle->getTime());
                    
                    h_HcalHit_Time_Matched_Tdif_SD[9]->Fill(part_hcalhit_timeDiff);
                    h_HcalHit_Time_Matched_Tdif_SD[ecal_sumESD]->Fill(part_hcalhit_timeDiff);

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
                        h_HcalHit_Time_Matched_Nucleons_SD[9]->Fill(matchedParticle->getTime());
//                        h_HCalhit_nucleon_time_vs_energy_SD[9]->Fill(
//                                matchedParticle->getTime(),
//                                matchedParticle->getEnergy());
                        h_HcalHit_Time_Matched_Nucleons_SD[ecal_sumESD]->Fill(matchedParticle->getTime());
//                        h_HCalhit_nucleon_time_vs_energy_SD[ecal_sumESD]->Fill(
//                                matchedParticle->getTime(),
//                                matchedParticle->getEnergy());
                    }
                    
                    h_Particle_Energy_Matched_SD[9]->Fill(matchedParticle->getEnergy());
                    h_Particle_Energy_Matched_SD[ecal_sumESD]->Fill(matchedParticle->getEnergy());
        
                    switch(pdgID) {
                        case 11:
                            h_HcalHit_ZbyR_Matched_Electron_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            h_HcalHit_ZbyR_Matched_Electron_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            break;
                        case 22:
                            h_HcalHit_ZbyR_Matched_Photon_SD[9]->Fill(hcalhit->getZ(), hcalhit_radialdist); 
                            h_HcalHit_ZbyR_Matched_Photon_SD[ecal_sumESD]->Fill(hcalhit->getZ(), hcalhit_radialdist);
//                            h_HCalhit_photon_energy_SD[9]->Fill(matchedParticle->getEnergy());
//                            h_HCalhit_photon_energy_SD[ecal_sumESD]->Fill(matchedParticle->getEnergy());
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
        
        numNonNoiseHits_ = 0;
        numMatchedHits_ = 0;
        numEvents_ = 0;

        // Make directory tree to organize histograms in
        TDirectory* d_base = getHistoDirectory();

        //////////// Directory Tree
        //  Ecal_SummedEnergy
        //  NumParticles
        //  EventMaxPE
        //  Particle
        //      PDGID
        //          All
        //          Matched
        //      HitDistance
        //          All
        //          Matched
        //      Energy
        //          All
        //          Matched
        //  HcalHit
        //      Z
        //      ZbyR
        //          All
        //          Unmatched
        //          TimeLess15
        //          TimeGreat40
        //          Matched
        //              Photon
        //              Electron
        //              Neutron
        //              Other
        //              TdifLess15
        //              TdifGreat40
        //      PE
        //          All
        //          TimeLess15
        //          TimeGreat40
        //          Matched
        //              TdifLess15
        //              TdifGreat40
        //      Time
        //          All
        //          Matched
        //              All
        //              Nucleons
        //              Tdif

        TDirectory* d_Ecal_SummedEnergy                 = d_base->mkdir( "Ecal_SummedEnergy" );
        TDirectory* d_NumParticles                      = d_base->mkdir( "NumParticles" );
        TDirectory* d_EventMaxPE                        = d_base->mkdir( "EventMaxPE" );

        TDirectory* d_Particle                          = d_base->mkdir( "Particle" );
        TDirectory* d_Particle_PDGID                    = d_Particle->mkdir( "PDGID" );
        TDirectory* d_Particle_PDGID_All                = d_Particle_PDGID->mkdir( "All" );
        TDirectory* d_Particle_PDGID_Matched            = d_Particle_PDGID->mkdir( "Matched" );
        TDirectory* d_Particle_HitDistance              = d_Particle->mkdir( "HitDistance" );
        TDirectory* d_Particle_HitDistance_All          = d_Particle_HitDistance->mkdir( "All" );
        TDirectory* d_Particle_HitDistance_Matched      = d_Particle_HitDistance->mkdir( "Matched" );
        TDirectory* d_Particle_Energy                   = d_Particle->mkdir( "Energy" );
        TDirectory* d_Particle_Energy_All               = d_Particle_Energy->mkdir( "All" );
        TDirectory* d_Particle_Energy_Matched           = d_Particle_Energy->mkdir( "Matched" );

        TDirectory* d_HcalHit                           = d_base->mkdir( "HcalHit" );
        TDirectory* d_HcalHit_Z                         = d_HcalHit->mkdir( "Z" );
        TDirectory* d_HcalHit_ZbyR                      = d_HcalHit->mkdir( "ZbyR" );
        TDirectory* d_HcalHit_ZbyR_All                  = d_HcalHit_ZbyR->mkdir( "All" );
        TDirectory* d_HcalHit_ZbyR_Unmatched            = d_HcalHit_ZbyR->mkdir( "Unmatched" );
        TDirectory* d_HcalHit_ZbyR_TimeLess15           = d_HcalHit_ZbyR->mkdir( "TimeLess15" );
        TDirectory* d_HcalHit_ZbyR_TimeGreat40          = d_HcalHit_ZbyR->mkdir( "TimeGreat40" );
        TDirectory* d_HcalHit_ZbyR_Matched              = d_HcalHit_ZbyR->mkdir( "Matched" );
        TDirectory* d_HcalHit_ZbyR_Matched_Photon       = d_HcalHit_ZbyR_Matched->mkdir( "Photon" );
        TDirectory* d_HcalHit_ZbyR_Matched_Electron     = d_HcalHit_ZbyR_Matched->mkdir( "Electron" );
        TDirectory* d_HcalHit_ZbyR_Matched_Neutron      = d_HcalHit_ZbyR_Matched->mkdir( "Neutron" );
        TDirectory* d_HcalHit_ZbyR_Matched_Other        = d_HcalHit_ZbyR_Matched->mkdir( "Other" );
        TDirectory* d_HcalHit_ZbyR_Matched_TdifLess15   = d_HcalHit_ZbyR_Matched->mkdir( "TdifLess15" );
        TDirectory* d_HcalHit_ZbyR_Matched_TdifGreat40  = d_HcalHit_ZbyR_Matched->mkdir( "TdifGreat40" );
        TDirectory* d_HcalHit_PE                        = d_HcalHit->mkdir( "PE" );
        TDirectory* d_HcalHit_PE_All                    = d_HcalHit_PE->mkdir( "All" );
        TDirectory* d_HcalHit_PE_TimeLess15             = d_HcalHit_PE->mkdir( "TimeLess15" );
        TDirectory* d_HcalHit_PE_TimeGreat40            = d_HcalHit_PE->mkdir( "TimeGreat40" );
        TDirectory* d_HcalHit_PE_Matched                = d_HcalHit_PE->mkdir( "Matched" );
        TDirectory* d_HcalHit_PE_Matched_TdifLess15     = d_HcalHit_PE_Matched->mkdir( "TdifLess15" );
        TDirectory* d_HcalHit_PE_Matched_TdifGreat40    = d_HcalHit_PE_Matched->mkdir( "TdifGreat40" );
        TDirectory* d_HcalHit_Time                      = d_HcalHit->mkdir( "Time" );
        TDirectory* d_HcalHit_Time_All                  = d_HcalHit_Time->mkdir( "All" );
        TDirectory* d_HcalHit_Time_Matched              = d_HcalHit_Time->mkdir( "Matched" );
        TDirectory* d_HcalHit_Time_Matched_All          = d_HcalHit_Time_Matched->mkdir( "All" );
        TDirectory* d_HcalHit_Time_Matched_Nucleons     = d_HcalHit_Time_Matched->mkdir( "Nucleons" );
        TDirectory* d_HcalHit_Time_Matched_Tdif         = d_HcalHit_Time_Matched->mkdir( "Tdif" );

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
        range_names[8] = "ninf-n6";
        range_names[9] = "ninf-pinf";

        for ( int i = 0; i < range_names.size(); i++ ) {
            std::string range_name = range_names.at(i);

            d_Ecal_SummedEnergy->cd();
            h_Ecal_SummedEnergy_SD[i]=new TH1D(
                    ("Ecal_SummedEnergy_SD_"+range_name).c_str(),
                    ("Ecal Summed Energy SD "+range_name+";Energy(MeV) (10MeV bin width);Count").c_str(),
                    800,0,8000);//10MeV bins

            d_NumParticles->cd();
            h_NumParticles_SD[i]=new TH1D(
                    ("NumParticles_SD_"+range_name).c_str(),
                    ("Num Particles in Event SD "+range_name+";Number of particles per event;Count").c_str(),
                    50,0,50);

            d_EventMaxPE->cd();
            h_EventMaxPE_SD[i]=new TH1D(
                    ("EventMaxPE_SD_"+range_name).c_str(),
                    ("Maximum PE for all Hits in Event SD "+range_name+";Photoelectrons(PEs);Count").c_str(),
                    500,0,500);

            d_Particle_PDGID_All->cd();
            h_Particle_PDGID_All_SD[i]=new TH1D(
                    ("Particle_PDGID_All_SD_"+range_name).c_str(),
                    ("PDG IDs SD "+range_name).c_str(),
                    10000,-5000,5000);

            d_Particle_PDGID_Matched->cd();
            h_Particle_PDGID_Matched_SD[i]=new TH1D(
                    ("Particle_PDGID_Matched_SD_"+range_name).c_str(),
                    ("PDG IDs when matched SD "+range_name).c_str(),
                    10000,-5000,5000);

            d_Particle_HitDistance_All->cd();
            h_Particle_HitDistance_All_SD[i]=new TH1D(
                    ("Particle_HitDistance_All_SD_"+range_name).c_str(), 
                    ("Distance between SimParticle and HcalHit SD "+range_name+" (5mm bins)").c_str(), 
                    400, 0, 2000);
            
            d_Particle_HitDistance_Matched->cd();
            h_Particle_HitDistance_Matched_SD[i]=new TH1D(
                    ("Particle_HitDistance_Matched_SD_"+range_name).c_str(), 
                    ("Distance between SimParticle and HcalHit when matched SD "+range_name+" (5mm bins)").c_str(), 
                    400, 0, 2000);
            
            d_Particle_Energy_All->cd();
            h_Particle_Energy_All_SD[i]=new TH1D(
                    ("Particle_Energy_All_SD_"+range_name).c_str(),
                    ("All Particle Energies SD "+range_name+";Energy(MeV)(5MeV bin width);Count").c_str(),
                    800,0,4000);

            d_Particle_Energy_Matched->cd();
            h_Particle_Energy_Matched_SD[i]=new TH1D(
                    ("Particle_Energy_Matched_SD_"+range_name).c_str(),
                    ("Matched Particle Energies SD "+range_name+";Energy(MeV)(5MeV bin width);Count").c_str(),
                    800,0,4000);

            d_HcalHit_Z->cd();
            h_HcalHit_Z_SD[i]=new TH1D(
                    ("HcalHit_Z_SD_"+range_name).c_str(),
                    ("Z depth of HCal hits SD "+range_name+" (10mm bins)").c_str(),
                    320, 0, 3200);

            d_HcalHit_ZbyR_All->cd();
            h_HcalHit_ZbyR_All_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_All_SD_"+range_name).c_str(), 
                    ("All HcalHit locations SD "+range_name+";Z depth (mm); radial distance from z-axis (mm)").c_str(),
                    80,0,3200,112,0,4500);

            d_HcalHit_ZbyR_Unmatched->cd();
            h_HcalHit_ZbyR_Unmatched_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Unmatched_SD_"+range_name).c_str(), 
                    ("Hcal unmatched hit locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            d_HcalHit_ZbyR_TimeLess15->cd();
            h_HcalHit_ZbyR_TimeLess15_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_TimeLess15_SD_"+range_name).c_str(),
                    ("HcalHits with Time < 15ns locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            d_HcalHit_ZbyR_TimeGreat40->cd();
            h_HcalHit_ZbyR_TimeGreat40_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_TimeGreat40_SD_"+range_name).c_str(),
                    ("HcalHits with Time > 40ns locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            d_HcalHit_ZbyR_Matched_Photon->cd();
            h_HcalHit_ZbyR_Matched_Photon_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Matched_Photon_SD_"+range_name).c_str(), 
                    ("Hcal photon hit locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            d_HcalHit_ZbyR_Matched_Electron->cd();
            h_HcalHit_ZbyR_Matched_Electron_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Matched_Electron_SD_"+range_name).c_str(), 
                    ("Hcal electron hit locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            d_HcalHit_ZbyR_Matched_Neutron->cd();
            h_HcalHit_ZbyR_Matched_Neutron_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Matched_Neutron_SD_"+range_name).c_str(), 
                    ("Hcal neutron hit locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            d_HcalHit_ZbyR_Matched_Other->cd();
            h_HcalHit_ZbyR_Matched_Other_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Matched_Other_SD_"+range_name).c_str(), 
                    ("Hcal other particle hit locations SD "+range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            d_HcalHit_ZbyR_Matched_TdifLess15->cd();
            h_HcalHit_ZbyR_Matched_TdifLess15_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Matched_TdifLess15_SD_"+range_name).c_str(),
                    ("Matched HcalHit location with time dif < 15ns SD "+
                     range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            d_HcalHit_ZbyR_Matched_TdifGreat40->cd();
            h_HcalHit_ZbyR_Matched_TdifGreat40_SD[i]=new TH2D(
                    ("HcalHit_ZbyR_Matched_TdifGreat40_SD_"+range_name).c_str(),
                    ("Matched HcalHit location with time dif > 40ns SD "+
                     range_name+";Z depth(mm);radial distance from z-axis(mm)").c_str(),
                    80,0,3200,112,0,4500);

            d_HcalHit_PE_All->cd();
            h_HcalHit_PE_All_SD[i]=new TH1D(
                    ("HcalHit_PE_All_SD_"+range_name).c_str(),
                    ("PEs of all HcalHits SD "+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            d_HcalHit_PE_TimeLess15->cd();
            h_HcalHit_PE_TimeLess15_SD[i]=new TH1D(
                    ("HcalHit_PE_TimeLess15_SD_"+range_name).c_str(),
                   ("HcalHits with Time < 15ns SD "+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            d_HcalHit_PE_TimeGreat40->cd();
            h_HcalHit_PE_TimeGreat40_SD[i]=new TH1D(
                    ("HcalHit_PE_TimeGreat40_SD_"+range_name).c_str(),
                    ("HcalHits with Time > 40ns SD "+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            d_HcalHit_PE_Matched_TdifLess15->cd();
            h_HcalHit_PE_Matched_TdifLess15_SD[i]=new TH1D(
                    ("HcalHit_PE_Matched_TdifLess15_SD_"+range_name).c_str(),
                    ("Matched HcalHit with Time dif < 15ns SD "+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            d_HcalHit_PE_Matched_TdifGreat40->cd();
            h_HcalHit_PE_Matched_TdifGreat40_SD[i]=new TH1D(
                    ("HcalHit_PE_Matched_TdifGreat40_SD_"+range_name).c_str(),
                    ("Matched HcalHit with Time dif > 40ns SD "+range_name+";Photoelectrons(PEs);Count").c_str(),
                    200,0,200);

            d_HcalHit_Time_All->cd();
            h_HcalHit_Time_All_SD[i]=new TH1D(
                    ("HcalHit_Time_All_SD_"+range_name).c_str(),
                    ("Time of All HcalHits SD "+range_name+";time(ns)(5ns bin width);Count").c_str(),
                    100,0,500);

            d_HcalHit_Time_Matched_All->cd();
            h_HcalHit_Time_Matched_All_SD[i]=new TH1D(
                    ("HcalHit_Time_Matched_All_SD_"+range_name).c_str(), 
                    ("Time of Matched HcalHits SD "+range_name+";Time(ns);Number of particles created").c_str(), 
                    500, 0, 500);

            d_HcalHit_Time_Matched_Nucleons->cd();
            h_HcalHit_Time_Matched_Nucleons_SD[i]=new TH1D(
                    ("HcalHit_Time_Matched_Nucleons_SD_"+range_name).c_str(),
                    ("Time of HcalHits Matched to Nucleons SD "+range_name+";Time(ns);Number of Nucleons created").c_str(), 
                    500, 0, 500);
            
            d_HcalHit_Time_Matched_Tdif->cd();
            h_HcalHit_Time_Matched_Tdif_SD[i]=new TH1D(
                    ("HcalHit_Time_Matched_Tdif_SD_"+range_name).c_str(),
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
        
        double hitRate, matchRate;
        {
            //temporary variables for calculating rates
            double numerator = numNonNoiseHits_;
            double denominator = numEvents_;
            hitRate = numerator/denominator;

            numerator = numMatchedHits_;
            denominator = numNonNoiseHits_;
            matchRate = numerator / denominator;
        }

        printf( "Number of Events:          %i\n" , numEvents_ );
        printf( "Number of Non Noise Hits:  %i\n" , numNonNoiseHits_ );
        printf( "Number of Matched Hits:    %i\n" , numMatchedHits_ );
        printf( "Hit Rate (hits/events):    %f\n" , hitRate );
        printf( "Match Rate (matches/hits): %f\n" , matchRate );

        return;
    }

} //ldmx namespace

DECLARE_ANALYZER_NS(ldmx, HcalHitMatcher);
