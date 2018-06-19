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
        
        //Log hits by layer
        std::vector< std::vector<ldmx::HcalHit*> > hitlog;

	    for (size_t i=0; i<tca->GetEntriesFast(); i++) {
            
            //get current hit
            const ldmx::HcalHit* chit=(const ldmx::HcalHit*)(tca->At(i));
            float curr_PE = chit->getPE();

            //only process non-noise hits
            if ( curr_PE > minPE_ ) {
                
                h_includedhits->Fill( curr_PE );

                //Only process if hit is in the back Hcal
                if ( chit->getSection() == ldmx::HcalSection::BACK ) {

                } //only process if hit is in back hcal
                else {
                    nNonBack_++;
                }

            } //only process non-noise hits
	    
        } //loop through all entries in calorimeter hits for current event (i)
        
        

    }

    void HcalLayerAnalyzer::onProcessStart() {
        getHistoDirectory();
        nNonBack_ = 0;
        h_includedhits = new TH1F("h_includedhits","PE Distribution of included hits",500,0.5,500.5);
        //Declare a histogram per layer of back hcal
        /*
        for ( ) {
            TH1F* h_newhistogram = new TH1F("name","title",500,0.0,1.0); //declare new histogram
            h_hitsperlayer_.pushback( h_newhistogram ); //put new histogram into vector
        } //loop through all layers in back hcal
        */
    }

    void HcalLayerAnalyzer::onProcessEnd() {
        std::cout << "Number Hits NOT in Back Hcal: " << nNonBack_ << std::endl;
    }

}

DECLARE_ANALYZER_NS(ldmx, HcalLayerAnalyzer);
