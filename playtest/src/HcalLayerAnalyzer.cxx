/**
 * @file HcalLayerAnalyzer.cxx
 * @brief Implementation file for HcalLayerAnalyzer
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "playtest/HcalLayerAnalyzer.h"

namespace ldmx {

    void HcalLayerAnalyzer::configure(const ldmx::ParameterSet& ps) {
        caloCol_=ps.getString("caloHitCollection");
        minPE_ = ps.getFloat("minPE");
    }

    void HcalLayerAnalyzer::analyze(const ldmx::Event& event) {
        const TClonesArray* tca=event.getCollection(caloCol_);
        
        //Log hits by layer
        std::vector< std::vector<ldmx::HcalHit*> > hitlog;

	    for (size_t i=0; i<tca->GetEntriesFast(); i++) {
            
            //get current hit
            const ldmx::HcalHit* chit=(const ldmx::HcalHit*)(tca->At(i));
            
            //only process non-noise hits
            if ( chit->getPE() > minPE_ ) {
                
                h_includedhits->Fill( chit->getPE() );

                //Only process if hit is in the back Hcal
                if ( chit->getSection() == ldmx::HcalSection::BACK ) {
                    
                    //obtain hit information
                    float curr_layer = chit->getLayer();
                    float curr_strip = chit->getStrip();
               
                    //Start iterator
                    std::vector< std::vector<ldmx::HcalHit*> >::iterator it = hitlog.begin();
                    while ( it != hitlog.end() ) {
                        if (curr_layer = it[0]->getLayer() ) {
                            it->push_back( chit );
                            break;
                        }

                        ++it;
                    }   

                    if ( it == hitlog.end() ) {
                        std::vector<ldmx::HcalHit*> tmp;
                        tmp.push_back(chit);
                        hitlog.push_back(tmp);
                    }

                } //only process if hit is in back hcal
                else {
                    nNonBack_++;
                }

            } //only process non-noise hits
	    
        } //loop through all entries in calorimeter hits for current event (i)
        
        //Sort hitlog
        std::sort( hitlog.begin() , hitlog.end() , compareLayer );
        for (std::vector< std::vector<ldmx::HcalHit*> >::iterator it = hitlog.begin(); it != hitlog.end(); ++it ) {
            std::sort( it->begin() , it->end() , compareStrip );
        } //loop through layers in hitlog (i)

        

    }

    void HcalLayerAnalyzer::onProcessStart() {
        getHistoDirectory();
        maxlayer_ = 0.0;
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

    bool HcalLayerAnalyzer::compareStrip( const ldmx::HcalHit* lhs , const ldmx::HcalHit* rhs ) {
        return ( lhs->getStrip() < rhs->getStrip() );
    }

    bool HcalLayerAnalyzer::compareLayer( const std::vector<ldmx::HcalHit*> lhs , const std::vector<ldmx::HcalHit*> rhs ) {
        return ( lhs[0]->getLayer() < rhs[0]->getLayer() );
    }
}

DECLARE_ANALYZER_NS(ldmx, HcalLayerAnalyzer);
