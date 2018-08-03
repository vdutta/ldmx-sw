/**
 * @file HcalMipTrackProducer.cxx
 * @brief Implementation file for HcalMipTrackProducer class.
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Hcal/HcalMipTrackProducer.h"

namespace ldmx {

    void HcalMipTrackProducer::configure(const ldmx::ParameterSet& ps) {

        hcalHitCollName_ = ps.getString( "HcalHitCollectionName" );

        hcalHitPassName_ = ps.getString( "HcalHitPassName" );

        TClonesArray *hcalMipTracks_ = new TClonesArray( EventConstants::HCAL_MIP_TRACK.c_str() , 10000 );

        hcalMipTracksCollName_ = ps.getString( "HcalMipTrackCollectionName" );

        trackRadius_ = ps.getDouble( "MipTrackRadius" );

        minPE_ = ps.getDouble( "MinimumPE" );

        maxEnergy_ = ps.getDouble( "MaximumEnergy" );

        minNumClusters_ = ps.getInteger( "MinimumNumClusters" );
        
        //tracks are started with a pair of points, so
        //  the absolute minimum number of clusters in a track is 2
        if ( minNumClusters_ < 2 ) {
            minNumClusters_ = 2;
        } //minNumClusters_ validation check
        
        return;
    }

    void HcalMipTrackProducer::produce(ldmx::Event& event) {
        
        //counts clusters and enumerates them
        int clustercnt = 0;

        const TClonesArray *rawhits = event.getCollection( hcalHitCollName_ , hcalHitPassName_ );

        //go through raw hits and throw away noise hits
        
        //cluster hits
        
        //iterate through all pairs of points
        //  construct track in cylinder
        //  check if plausible track
        //  create HcalMipTrack instance

        //get best HcalMipTrack

        //store best in collection
        //delete best from cluster log

        //repeat track construction until no more pairs of clusters

        return;
    }
    
}

DECLARE_PRODUCER_NS(ldmx, HcalMipTrackProducer);
