/**
 * @file TriggerAnalyzer.cxx
 * @brief Implementation file for TriggerAnalyzer class
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/TriggerAnalyzer.h"

namespace ldmx {

    void TriggerAnalyzer::configure(const ldmx::ParameterSet& ps) {
        
        hcalTriggerObjectName_ = ps.getString( "HcalTriggerObjectName" );

        hcalTriggerPassName_ = ps.getString( "HcalTriggerPassName" );

        numFalsePass_ = 0;
        numTruePass_ = 0;
        numFalseFail_ = 0;
        numTrueFail_ = 0;

        return;
    }

    void TriggerAnalyzer::analyze(const ldmx::Event& event) {
        
        //get list of triggers
        const TClonesArray *triggers = event.getCollection( "Trigger" , hcalTriggerPassName_ );
        //get list of actual sim particles
        const TClonesArray *simparticles = event.getCollection( "SimParticles" , "sim" );
        //get number of hcalHits that aren't noise
        const TClonesArray *hcalHits = event.getCollection( "hcalDigis" , "recon" );
        int nHcalHits = 0;
        for ( int iH = 0; iH < hcalHits->GetEntriesFast(); iH++ ) {
            HcalHit *chit = (HcalHit *)(hcalHits->At(iH));
            if ( !chit->getNoise() )
                nHcalHits++;
        }
        
        //get hcal  trigger
        int ntriggers = triggers->GetEntriesFast();
        const TriggerResult *hcalTrigger;
        for ( int i = 0; i < ntriggers; i++ ) {
            hcalTrigger = (const TriggerResult *)(triggers->At(i));
            if ( hcalTrigger->getName() == hcalTriggerObjectName_ )
                break;
        }

        //count number of actual muons
        int nmuons = 0;
        for ( int iP = 0; iP < simparticles->GetEntriesFast(); iP++ ) {
            SimParticle *simp = (SimParticle *)(simparticles->At(iP));
            if ( simp->getPdgID() == 13 or simp->getPdgID() == -13 ) {
                nmuons++;
            } //if a muon
        }
        
        //ignore events with zero hcalHits (boring)
        if ( nHcalHits > 0 ) {
            //check accuracy of  trigger and set storage hint
            bool triggerpass = hcalTrigger->passed();
            bool realpass = ( nmuons > 0 );
            if ( triggerpass and realpass) {
                numTruePass_++;
                setStorageHint( hint_mustDrop );
            } else if ( triggerpass and !realpass ) {
                numFalsePass_++;
                setStorageHint( hint_mustKeep );
            } else if ( !triggerpass and realpass ) {
                numFalseFail_++;
                setStorageHint( hint_mustKeep );
            } else {
                numTrueFail_++;
                setStorageHint( hint_mustDrop );
            }
        } else {
            setStorageHint( hint_mustDrop );
        }
        
        return;
    }

    void TriggerAnalyzer::onProcessEnd() {
        
        unsigned int numEvents = numTruePass_+numTrueFail_+numFalsePass_+numFalseFail_;
        double accuracy = ((double)(numTruePass_) + (double)(numTrueFail_))/(double)(numEvents);
        double sensitivity = (double)(numTruePass_)/(double)(numTruePass_ + numFalseFail_);
        double precision = (double)(numTruePass_)/(double)(numTruePass_ + numFalsePass_);
        double missrate = 1 - sensitivity;
        double falsePassRate = (double)(numFalsePass_)/(double)(numFalsePass_+numTruePass_);
        printf( "\n" );
        printf( " ===============================\n" );
        printf( " | %27s |\n" , hcalTriggerObjectName_.c_str() );
        printf( " |       Confusion Table       |\n" );
        printf( " |         ||    Sim Particle  |\n" );
        printf( " | Trigger ||   Pass | Fail    |\n" );
        printf( " |    Pass ||%7d | %-7d |\n" , numTruePass_ , numFalsePass_ );
        printf( " |    Fail ||%7d | %-7d |\n" , numFalseFail_ , numTrueFail_ );
        printf( " |=============================|\n" );
        printf( " | N Events        | %-9d |\n" , numEvents );
        printf( " | Accuracy        | %-9f |\n" , accuracy );
        printf( " | True Pass Rate  | %-9f |\n" , sensitivity );
        printf( " | False Fail Rate | %-9f |\n" , missrate );
        printf( " | False Pass Rate | %-9f |\n" , falsePassRate );
        printf( " | Precision       | %-9f |\n" , precision );
        printf( " | Informedness    | %-9f |\n" , sensitivity - falsePassRate );
        printf( " ===============================\n" );

        return;
    }

}

DECLARE_ANALYZER_NS(ldmx, TriggerAnalyzer);
