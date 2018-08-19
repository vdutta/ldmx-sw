/**
 * @file MuonAnalyzer.cxx
 * @brief Implementation file for MuonAnalyzer class.
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/MuonAnalyzer.h"

namespace ldmx {

    void MuonAnalyzer::configure(const ldmx::ParameterSet& ps) {

        return;
    }

    void MuonAnalyzer::analyze(const ldmx::Event& event) {
        
        //one set for each orientation
        std::set<int> layersHit[5], stripsHit[5];
        
        //get hcal hits
        const TClonesArray *hcalHits = event.getCollection( "hcalDigis" , "recon" );
        
        //fill sets of layers and strips (sorting them)
        int nhits = hcalHits->GetEntriesFast();
        for ( int iH = 0; iH < nhits; iH++ ) {
            const HcalHit *hit = (const HcalHit *)(hcalHits->At(iH));
            int section = hit->getSection();
            stripsHit[section].insert( hit->getStrip() );
            layersHit[section].insert( hit->getLayer() );
        }

        //fill histograms
        for ( int s = 0; s < 5; s++ ) {
            if ( !stripsHit[s].empty() ) {
                hMinStrip_[s]->Fill( *stripsHit[s].begin() );
                hMaxStrip_[s]->Fill( *std::prev(stripsHit[s].end()) );
                hNumConsecStrips_[s]->Fill( consecutiveCount( stripsHit[s] ) );
            }
    
            if ( !layersHit[s].empty() ) {
                hMinLayer_[s]->Fill( *layersHit[s].begin() );
                hMaxLayer_[s]->Fill( *std::prev(layersHit[s].end()) );
                hNumConsecLayers_[s]->Fill( consecutiveCount( layersHit[s] ) );
            }
        }

        return;
    }

    void MuonAnalyzer::onProcessStart() {
        
        getHistoDirectory();
        
        for ( int s = 0; s < 5; s++ ) {

            std::string sect = std::to_string(s);

            hMinStrip_[s] = new TH1F( ("hMinStrip_"+sect).c_str() , ("Minimum Strips in Hcal Section "+sect).c_str() , 
                40 , 0.0 , 40.0 );
            hMaxStrip_[s] = new TH1F( ("hMaxStrip_"+sect).c_str() , ("Maximum Strips in Hcal Section "+sect).c_str() ,
                40 , 0.0 , 40.0 );
            hNumConsecStrips_[s] = new TH1F( ("hNumConsecStrips_"+sect).c_str() , ("Num Consecutive Strips in Hcal Section "+sect).c_str() ,
                40 , 0.0 , 40.0 );
            
            hMinLayer_[s] = new TH1F( ("hMinLayer_"+sect).c_str(), ("Minimum Layers in Hcal Section "+sect).c_str(),
                150 , 0.0 , 150.0 );
            hMaxLayer_[s] = new TH1F( ("hMaxLayer_" +sect).c_str() , ("Maximum Layers in Hcal Section "+sect).c_str() ,
                150 , 0.0 , 150.0 );
            hNumConsecLayers_[s] = new TH1F( ("hNumConsecLayers_"+sect).c_str() , ("Num Consecutive Layers in Hcal Section "+sect).c_str() ,
                150 , 0.0 , 150.0 );
        }
        
        
        return;
    }
    
    int MuonAnalyzer::consecutiveCount( const std::set<int> &list ) const {
        
        int consec(0), maxConsec(0);
        
        if ( !list.empty() ) {
            
            int prev = *list.cbegin();
            for ( int curr : list ) {
                   
                if ( curr - prev > 1 ) {
                    if ( consec > maxConsec ) //improve upon max
                        maxConsec = consec;
    
                    consec = 0; //reset counter
                } //check if not consecutive anymore
     
                consec++; //count this current number

                prev = curr; //update prev
            } //go through ordered list
    
            if ( consec > maxConsec )
                maxConsec = consec;
    
        } //list is non-empty

        return maxConsec;

    }
}

DECLARE_ANALYZER_NS(ldmx, MuonAnalyzer);
