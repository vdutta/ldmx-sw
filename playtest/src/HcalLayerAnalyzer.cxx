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
        nStrips_= ps.getInteger("nStrips");
        nEcalThickness_ = ps.getInteger("nEcalThickness");

        origin_ = static_cast<float>(nStrips_)/2;
        lowside_ = origin_ - static_cast<float>(nEcalThickness_)/2;
        upside_ = origin_ + static_cast<float>(nEcalThickness_)/2;
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

        for (HitLog::iterator it = log.begin(); it != log.end(); ++it) {
            std::cout << it->first << std::endl; //printing out list of keys
        }

    }

    void HcalLayerAnalyzer::onProcessStart() {
        getHistoDirectory();
        h_includedhits = new TH1F("h_includedhits","PE Distribution of included hits",500,0.5,500.5);

        nNotIncluded_ = 0;
        layermod_ = 1000;
    }

    void HcalLayerAnalyzer::onProcessEnd() {
        std::cout << "Number Hits NOT included in analysis: " << nNotIncluded_ << std::endl;
    }

    int HcalLayerAnalyzer::keygen( HitPtr hit ) const {
        return static_cast<int>( hit->getLayer()*layermod_ + hit->getStrip() );
    }
    

    std::pair< HitPtr , HitPtr > HcalLayerAnalyzer::search( const HitLog log , const int lowkey , const int upkey ) const {
        
        std::pair< HitPtr , HitPtr > ret( nullptr , nullptr );
        
        if ( lowkey > upkey ) { //Mis-use correction
            std::cout << "Input strip numbers to search in wrong order: " << lowstrip << " " << upstrip << std::endl;
            std::cout << "Returning an empty search" << std::endl;
            return ret;
        } else { //inputs are correct form

            auto lowbound = log.lower_bound( lowkey ); //points to first key that is not before lowkey (equivalent or after) (map::end if all are before lowkey)
            auto upbound = log.upper_bound( upkey ); //points to first key after upkey (map::end if nothing after upkey)
        
            if ( lowbound != upbound ) { //check to see if range has any thickness
                //there is at least one hit in the range
                //lowbound points to the first one
                //see if there is a hit on either side of lowerbound
                auto beforeside = std::prev( lowbound ); 
                auto afterside = std::next( lowbound );
                
                int beforekeydif = lowbound->first - beforeside->first;
                int afterkeydif = afterside->first - lowbound->first;

                if ( beforekeydif != 1 or afterkeydif != 1 ) {
                    //lowbound has at most one neighbor
                    ret.first = lowbound->second;
                    
                    if ( beforekeydif == 1 ) {
                        //beforeside is the neighbor for lowbound
                        ret.second = beforeside->second;
                    } else if ( afterkeydif == 1 ) {
                        //afterside is the neighbor for lowerbound
                        ret.second = afterside->second;
                    } //else: lowbound is truly isolated

                } //check if lowbound could be isolated

            } //check to see if range has any thickness

        } //make sure inputs are correct
            
        return ret;
    }

    bool HcalLayerAnalyzer::stripbounds( const int seedlayer , const int seedstrip , const int layer , int &lowstrip , int &upstrip ) const {
        
        float slope = (seedstrip - origin_)/seedlayer; //slope of line between seed and origin

        lowstrip = static_cast<int>(floor( layer*slope + lowside_ ));
        upstrip = static_cast<int>(ceil( layer*slope + upside_ ));
        
        //Change if out of bounds
        if ( lowstrip < 1 )
            lowstrip = 1;
        if ( upstrip > nStrips_ )
            upstrip = nStrips_;

        if ( lowstrip > nStrips_ or upstrip < 1 ) {
            return false; //projected track is outside of Hcal
        }

        return true;
    }

    bool HcalLayerAnalyzer::findseed( const HitLog log , int &seedlayer , int &seedstrip ) const {
        
        std::pair< HitPtr , HitPtr > result = search( log , seedlayer*layermod_ , (seedlayer+1)*layermod_ - 1 );

        if ( result.first != nullptr ) { //Check if search was successful
            seedstrip = result.first->getStrip();
        } else {
            std::cout << "Warning:[HcalLayerAnalzyer::findseed] Unable to find isolated hit in input seedlayer." << std::endl;
            std::cout << "\tProceeding with first isolated hit found." << std::endl;
            result = search( log , 0 , 100*layermod_ );

            if ( result.first == nullptr ) {
                std::cout << "Error:[HcalLayerAnalyzer::findseed] Unable to find isolated hit anywhere in this event." << std::endl;
                return false;
            }

            seedlayer = result.first->getLayer();
            seedstrip = result.first->getStrip();
        }

        return true;
    }

}

DECLARE_ANALYZER_NS(ldmx, HcalLayerAnalyzer);
