#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "EventProc/TrigScintDigiProducer.h"

#include <iostream>
#include <exception>

namespace ldmx {

    TrigScintDigiProducer::TrigScintDigiProducer(const std::string& name, Process& process) :
        Producer(name, process) {
        hits_ = new TClonesArray("ldmx::TrigScintHit");
    }

    void TrigScintDigiProducer::configure(const ParameterSet& ps) {
        detID_       = new DefaultDetectorID();
        random_       = new TRandom3(ps.getInteger("randomSeed", 1000));
        NUM_STRIPS_PER_ARRAY_      = ps.getInteger("number_of_strips");
        NUM_ARRAYS_                = ps.getInteger("number_of_arrays");
        meanNoise_                 = ps.getDouble("meanNoise");
        mev_per_mip_               = ps.getDouble("mev_per_mip");
        pe_per_mip_                = ps.getDouble("pe_per_mip");
	input_collection_          = ps.getString("input_collection","TriggerPadUpSimHits");
	output_collection_          = ps.getString("output_collection","trigScintDigis");
    }

    unsigned int TrigScintDigiProducer::generateRandomID(TrigScintSection sec){
        DefaultDetectorID tempID;
        if( sec < TrigScintSection::NUM_SECTIONS ){
	    tempID.setFieldValue(0,int(sec));
	    tempID.setFieldValue(1,random_->Integer(NUM_STRIPS_PER_ARRAY_));
	}else
	    std::cout << "WARNING [TrigScintDigiProducer::generateRandomID]: TrigScintSection is not known" << std::endl;

        return tempID.pack();
    }

    void TrigScintDigiProducer::produce(Event& event) {

        std::map<unsigned int, int>   cellPEs;
        std::map<unsigned int, int>   cellMinPEs;
        std::map<unsigned int, float> Xpos, Ypos, Zpos, Edep, Time;
        std::unordered_set<unsigned int> noiseHitIDs;

        // looper over sim hits and aggregate energy depositions for each detID
	// AJW: what is the colletion type going to be?
        TClonesArray* simHits = (TClonesArray*) event.getCollection(input_collection_, "sim");

        int numSimHits = simHits->GetEntries();
        for (int iHit = 0; iHit < numSimHits; iHit++) {            
            SimCalorimeterHit* simHit = (SimCalorimeterHit*) simHits->At(iHit);
             int detIDraw = simHit->getID();
            detID_->setRawValue(detIDraw);
            detID_->unpack();
	    std::vector<float> position = simHit->getPosition();

            if (verbose_) {
                std::cout << "section: " << detID_->getFieldValue("section") << "  layer: " << detID_->getFieldValue("layer") <<  "  strip: " << detID_->getFieldValue("strip") <<std::endl;
            }        
            
            // for now, we take am energy weighted average of the hit in each stip to simulate the hit position. 
            // will use strip TOF and light yield between strips to estimate position.            
            if (Edep.find(detIDraw) == Edep.end()) {
                // first hit, initialize
                Edep[detIDraw] = simHit->getEdep();
                Time[detIDraw] = simHit->getTime() * simHit->getEdep();
                Xpos[detIDraw]      = position[0]* simHit->getEdep();
                Ypos[detIDraw]      = position[1]* simHit->getEdep();
                Zpos[detIDraw]      = position[2]* simHit->getEdep();
            } else {
                // not first hit, aggregate, and store the largest radius hit
                Xpos[detIDraw]      += position[0]* simHit->getEdep();
                Ypos[detIDraw]      += position[1]* simHit->getEdep();
                Zpos[detIDraw]      += position[2]* simHit->getEdep();
                Edep[detIDraw]      += simHit->getEdep();
		// AJW: need to figure out a better way to model this...
                Time[detIDraw]      += simHit->getTime() * simHit->getEdep();
            }
	    
        }

        // loop over detIDs and simulate number of PEs
        int ihit = 0;        
        for (std::map<unsigned int, float>::iterator it = Edep.begin(); it != Edep.end(); ++it) {
            int detIDraw = it->first;
            double depEnergy    = Edep[detIDraw];
            Time[detIDraw]      = Time[detIDraw] / Edep[detIDraw];
            Xpos[detIDraw]      = Xpos[detIDraw] / Edep[detIDraw];
            Ypos[detIDraw]      = Ypos[detIDraw] / Edep[detIDraw];
            Zpos[detIDraw]      = Zpos[detIDraw] / Edep[detIDraw];
            double meanPE       = depEnergy / mev_per_mip_ * pe_per_mip_;

	    //  dropping readout threshold for now... AJW
            //if( cellPEs[detIDraw] >= readoutThreshold_ ){ // > or >= ?
                
	  TrigScintHit *hit = (TrigScintHit*) (hits_->ConstructedAt(ihit));
          
	  hit->setID(detIDraw);
	  hit->setPE(cellPEs[detIDraw]);
	  hit->setMinPE(cellMinPEs[detIDraw]);
	  hit->setAmplitude(cellPEs[detIDraw]);
	  hit->setEnergy(depEnergy);
	  hit->setTime(Time[detIDraw]);
	  hit->setXpos(Xpos[detIDraw]); // quantized and smeared positions
	  hit->setYpos(Ypos[detIDraw]); // quantized and smeared positions
	  hit->setZpos(Zpos[detIDraw]);
	  hit->setNoise(false);
	  ihit++;
          
	  //}

            if (verbose_) {
                detID_->setRawValue(detIDraw);
                detID_->unpack();
                int layer = detID_->getFieldValue("layer");
                int subsection = detID_->getFieldValue("section");
                int strip = detID_->getFieldValue("strip");

                std::cout << "detID: " << detIDraw << std::endl;
                std::cout << "Layer: " << layer << std::endl;
                std::cout << "Subsection: " << subsection << std::endl;
                std::cout << "Strip: " << strip << std::endl;
                std::cout << "Edep: " << Edep[detIDraw] << std::endl;
                std::cout << "numPEs: " << cellPEs[detIDraw] << std::endl;
                std::cout << "time: " << Time[detIDraw] << std::endl;
                std::cout << "z: " << Zpos[detIDraw] << std::endl;
                std::cout << "Layer: " << layer << "\t Strip: " << strip << "\t X: " << Xpos[detIDraw] <<  "\t Y: " << Ypos[detIDraw] <<  "\t Z: " << Zpos[detIDraw] << std::endl;
            }        // end verbose            
        } 
        

        // ------------------------------- Noise simulation -------------------------------
	
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	
        event.add(output_collection_, hits_);
    }

}

DECLARE_PRODUCER_NS(ldmx, TrigScintDigiProducer);

