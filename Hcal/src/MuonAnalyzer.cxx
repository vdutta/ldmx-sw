/**
 * @file MuonAnalyzer.cxx
 * @brief Implementation file for MuonAnalyzer class.
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/MuonAnalyzer.h"

namespace ldmx {

    void MuonAnalyzer::configure(const ldmx::ParameterSet& ps) {
        
        triggerObjectName_ = ps.getString( "TriggerObjectName" );
        triggerPassName_ = ps.getString( "TriggerPassName" );
        
        hcalHitCollName_ = ps.getString( "HcalHitCollectionName" );
        hcalHitPassName_ = ps.getString( "HcalHitPassName" );

        return;
    }

    void MuonAnalyzer::analyze(const ldmx::Event& event) {
        
        //get trigger array
        const TClonesArray *triggers = event.getCollection( "Trigger" , triggerPassName_ );

        const TriggerResult *muonTrigger;
        int nTriggers = triggers->GetEntriesFast();
        int iT;
        for ( iT = 0; iT < nTriggers; iT++ ) {
            muonTrigger = (const TriggerResult *)(triggers->At(iT));
            if ( muonTrigger->getName() == triggerObjectName_ )
                break;
        }

        if ( iT >= nTriggers ) {
            std::cerr << "WARNING [ MuonAnalyzer ] : Could not find trigger object '" << triggerObjectName_ << "'." << std::endl;
        } else {
    
            //fill histograms for each section
            int maxconsecLayers(0), maxconsecStrips(0);
            for ( int s = 0; s < 5; s++ ) {
                int consecLayers = muonTrigger->getAlgoVar( 4*s+2 );
                int consecStrips = muonTrigger->getAlgoVar( 4*s+3 );
                
                hNumConsecLayers_[s]->Fill( consecLayers );
                hNumConsecStrips_[s]->Fill( consecStrips );
                
                if ( consecLayers > maxconsecLayers )
                    maxconsecLayers = consecLayers;

                if ( consecStrips > maxconsecStrips )
                    maxconsecStrips = consecStrips;
            }

        } //check if found trigger

        /*
         * Calculation of Path Length Uncertainty
         */

        //get hcal hits
        const TClonesArray *hcalHits = event.getCollection( hcalHitCollName_ , hcalHitPassName_ );

        //sort layers hit and strips hit into sections
        int nHits = hcalHits->GetEntriesFast();
        HcalHit *first(nullptr), *zero(nullptr);
        double dist2first(0);
        if ( nHits > 0 ) {
            zero = (HcalHit *)(hcalHits->At(0));
        }
        for ( int iH = 0; iH < nHits; iH++ ) {
            HcalHit *chit = (HcalHit *)(hcalHits->At(iH));
                
            if ( !chit->getNoise() and zero ) {
                
                double dx = zero->getX() - chit->getX();
                double dy = zero->getY() - chit->getY();
                double dz = zero->getZ() - chit->getZ();
                double dist = sqrt( dx*dx + dy*dy + dz*dz );

                if ( dist > dist2first ) {
                    dist2first = dist;
                    first = chit;
                }

            } //chit is not noise

        } //add each hit information to set

        //find last
        HcalHit *last(nullptr);
        if ( first ) {
            double dist2last(0);
            for ( int iH = 0; iH < nHits; iH++ ) {
                HcalHit *chit = (HcalHit *)(hcalHits->At(iH));
    
                if ( !chit->getNoise() ) {
    
                    double dx = first->getX() - chit->getX();
                    double dy = first->getY() - chit->getY();
                    double dz = first->getZ() - chit->getZ();
                    
                    double dist = sqrt( dx*dx + dy*dy + dz*dz );

                    if ( dist > dist2last ) {
                        dist2last = dist;
                        last = chit;
                    } //improve on distance
                } //ignore noise hits
            } //loop through all hits
        } //first exists
        
        double pathUnc(-1);
        if ( first != last ) {
            double dx = abs(last->getX() - first->getX());
            double dy = abs(last->getY() - first->getY());
            double dz = abs(last->getZ() - first->getZ());
            double dx2 = dx*dx;
            double dy2 = dy*dy;
            double dz2 = dz*dz;
            double scint_width = 100.0;
            double scint_thick = 6.0;
            double w2 = scint_width*scint_width;
            double t2 = scint_thick*scint_thick;
                
            if ( first->getSection() == 0 or last->getSection() == 0 ) {
                //back - layers along z
                double path = sqrt( t2*( 2 + dx2/dz2 + dy2/dz2 ) );
                pathUnc = sqrt( ( t2*(dx2+dy2)*(t2*(dx2+dy2)+w2*dz2) )/( sqrt(3)*dz2*dz2*(dx2+dy2+2*dz2) ) );
                pathUnc /= path;
            } else if ( first->getSection() < 3 and last->getSection() < 3 ) {
                //Top and/or bottom - layers along y
                double path = sqrt( t2*( 2 + dx2/dy2 + dz2/dy2 ) );
                pathUnc = sqrt( ( t2*(dx2+dz2)*(t2*(dx2+dz2)+w2*dy2) )/( sqrt(3)*dy2*dy2*(dx2+dz2+2*dy2) ) );
                pathUnc /= path;
            } else if ( first->getSection() > 2 and last->getSection() > 2 ) {
                //Left and/or right - layers along x
                double path = sqrt( t2*( 2 + dz2/dx2 + dy2/dx2 ) );
                pathUnc = sqrt( ( t2*(dz2+dy2)*(t2*(dz2+dy2)+w2*dx2) )/( sqrt(3)*dx2*dx2*(dz2+dy2+2*dx2) ) );
                pathUnc /= path;
            }
        }

        if ( pathUnc > 0.0 ) {
            hPathLengthUnc_->Fill( pathUnc );

            if ( muonTrigger->passed() )
                hPathLengthUncPassed_->Fill( pathUnc );

        }

        return;
    }

    void MuonAnalyzer::onProcessStart() {
        
        getHistoDirectory();
        
        for ( int s = 0; s < 5; s++ ) {

            std::string sect = std::to_string(s);

            hNumConsecStrips_[s] = new TH1F( ("hNumConsecStrips_"+triggerObjectName_+sect).c_str() , ("Num Consecutive Strips in Hcal Section "+sect).c_str() ,
                40 , 0.0 , 40.0 );
            hNumConsecLayers_[s] = new TH1F( ("hNumConsecLayers_"+triggerObjectName_+sect).c_str() , ("Num Consecutive Layers in Hcal Section "+sect).c_str() ,
                150 , 0.0 , 150.0 );
        }

        hPathLengthUnc_ = new TH1F( "hPathLengthUnc_" , "Uncertainty in Path Length for All Muons" ,
            100 , 0.0 , 1.0 );

        hPathLengthUncPassed_ = new TH1F( "hPathLengthUncPassed_" , "Uncertainty in Path Length for Passed Muons" ,
            100 , 0.0 , 1.0 );

        return;
    }
}

DECLARE_ANALYZER_NS(ldmx, MuonAnalyzer);
