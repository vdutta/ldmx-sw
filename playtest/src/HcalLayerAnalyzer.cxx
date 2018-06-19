/**
 * @file HcalLayerAnalyzer.cxx
 * @brief Implementation file for HcalLayerAnalyzer
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "playtest/HcalLayerAnalyzer.h"

namespace ldmx {
    
    void HcalLayerAnalyzer::configure(const ldmx::ParameterSet& ps) {
        caloCol_=ps.getString("caloHitCollection");
        minPE_ = static_cast<float>(ps.getDouble("minPE"));
    }

    void HcalLayerAnalyzer::analyze(const ldmx::Event& event) {
        const TClonesArray* tca=event.getCollection(caloCol_);
        
        HitLog log;

	    for (size_t i = 0; i < tca->GetEntriesFast(); i++) {
            
            //get current hit
            const ldmx::HcalHit* chit=(const ldmx::HcalHit*)(tca->At(i));
            float curr_PE = chit->getPE();

            //only process non-noise hits in the back hcal
            if ( curr_PE > minPE_ and chit->getSection() == ldmx::HcalSection::BACK ) {
                
                h_includedhits->Fill( curr_PE );

                int curr_key = keygen( chit );
                log[ curr_key ] = chit;

            } //only process non-noise hits in the back hcal
            else {
                nNotIncluded_++;
            }
	    
        } //loop through all entries in calorimeter hits for current event (i)
        
        //Now log has the non-noise hits in it
        //

        for (HitLog::iterator it = log.begin(); it != log.end(); ++it) {
            std::cout << it->first << std::endl; //printing out list of keys
        }

    }

    void HcalLayerAnalyzer::onProcessStart() {
        getHistoDirectory();
        h_includedhits = new TH1F("h_includedhits","PE Distribution of included hits",500,0.5,500.5);
        //Declare a histogram per layer of back hcal
        /*
        for ( ) {
            TH1F* h_newhistogram = new TH1F("name","title",500,0.0,1.0); //declare new histogram
            h_hitsperlayer_.pushback( h_newhistogram ); //put new histogram into vector
        } //loop through all layers in back hcal
        */

        nNotIncluded_ = 0;
        layermod_ = 1000;
    }

    void HcalLayerAnalyzer::onProcessEnd() {
        std::cout << "Number Hits NOT included in analysis: " << nNotIncluded_ << std::endl;
    }

    int HcalLayerAnalyzer::keygen( const ldmx::HcalHit* hit ) const {
        return static_cast<int>( hit->getLayer()*layermod_ + hit->getStrip() );
    }

}

DECLARE_ANALYZER_NS(ldmx, HcalLayerAnalyzer);
