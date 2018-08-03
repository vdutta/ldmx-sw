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
        
        const TClonesArray *rawhits = event.getCollection( hcalHitCollName_ , hcalHitPassName_ );

        //go through raw hits and ignore noise hits
        int nhits = rawhits.GetEntriesFast();
        for ( int iH = 0; iH < nhits; iH++ ) {
            
            HcalHit* chit = dynamic_cast<HcalHit*>((rawhits->At( iH ));

            if ( isNotNoise( chit ) ) {
                int section = chit->getSection();
                int layer = chit->getLayer();
                int strip = chit->getStrip();

                unsigned int key = section*1000*100 + layer*100 + strip;

                hcalHitLog_[ key ] = chit;
            } //if not noise hit

        } //iterate through rawhits (iH)

        clusterHits();
       
                //store best in collection
        //delete best from cluster log

        //repeat track construction until no more pairs of clusters

        return;
    }

    bool HcalMipTrackProducer::isNotNoise( const HcalHit* hit ) const {
        return ( !chit->getNoise() and chit->getPE() > minPE_ );
    }

    bool HcalMipTrackProducer::isMip( const MipCluster* cluster ) const {
        return ( cluster->getEnergy() < maxEnergy_ );
    }

    void HcalMipTrackProducer::clusterHits() {
        
        // cluster hits in same layer
        std::map< unsigned int , HcalHit* >::iterator itH;
        std::map< unsigned int , HCalHit* >::iterator prev_itH = hcalHitLog_.begin();
        MipCluster *current_cluster = new MipCluster();
        for ( itH = hcalHitLog_.begin(); itH != hcalHitLog_.end(); ++itH ) {
            
            //itH and prev_itH will both point to hcalHitLog_.begin() on first loop
            unsigned int keydif = itH->first - prev_itH->first;

            if ( keydif > 1 ) {
                //current hit is in different cluster
                // add previous cluster to log and clear temporary cluster
                current_cluster->setUID( prev_itH->first );
                current_cluster->set(); //calculate real space point
                if ( isMip( current_cluster ) ) {
                    clusterLog_[ current_cluster->getUID() ] = *current_cluster;
                }//check if mip
                delete current_cluster;
                current_cluster = new MipCluster();
            }//check if separate cluster

            current_cluster->addHit( itH->second );

            itH++;
            prev_itH = std::prev( itH );

        }//iterate through sorted hcalHitLog (itH)

        //clean up at end of hit log
        if ( current_cluster->getNumHits() > 0 ) {
            current_cluster->setUID( prev_itH->first );
            current_cluster->set();
            clusterLog_[ current_cluster->getUID() ] = *current_cluster;
        }

        if ( current_cluster ) {
            delete current_cluster;
            current_cluster = nullptr;
        }

        // cluster across layers (if overlap)
     
    }
    
    double HcalMipTrackProducer::distToLine( const std::vector<double> &P1 , const std::vector<double> &P2 ,
        const std::vector<double> &Q ) const {

        //calculate determinant of matrix with the three points as rows
        double det = P1[0]*P2[1]*Q[2] + P2[0]*Q[1]*P1[2] + Q[0]*P1[1]*P2[2] + P1[0]*P2[1]*Q[2]
            - Q[0]*P2[1]*P1[2] - P2[0]*P1[1]*Q[2] - P1[0]*Q[1]*P2[2] - Q[0]*P2[1]*P1[2];

        //distance between P1 and P2
        double dist = sqrt( (P2[0]-P1[0])*(P2[0]-P1[0]) + (P2[1]-P1[1])*(P2[1]-P1[1]) + (P2[2]-P1[2])*(P2[2]-P1[2]) );

        return ( std::abs(det)/dist );

    }

    bool HcalMipTrackProducer::compMipTracks( const HcalMipTrack &track1 , const HcalMipTrack &track2 ) const {
        
        bool better = false;
        if ( track1.isEmpty() ) {
            better = true;
        } else {
            //TEMPORARY
            better = ( track1.getEnergy() < track2.getEnergy() );
        } //if track1 is empty

        return better;
    }

    bool HcalMipTrackProducer::buildTrack( std::vector< unsigned int > track_mipids ) {
        //for clusters:
        //  no suffix means that it isn't an endpoint
        //  suffice {1,2} means that it is one of the endpoints
        bool success = false;
        track_mipids.clear();
        //iterate through all pairs of points
        HcalMipTrack best_track;
        std::map< unsigned int , MipCluster >::iterator itC1, itC2, itC; //iterators for map
        for ( itC1 = clusterLog_.begin(); itC1 != clusterLog.end(); ++itC1 ) {
            for ( itC2 = itC1+1; itC2 != clusterLog_.end(); ++itC2 ) {
                //construct track in cylinder
                
                std::vector< unsigned int > ctrack_mipids;

                std::vector< double > point1, point2 , errors1 , errors2;
                (itC1->second).getPoint( point1 , errors1 );
                (itC2->second).getPoint( point2 , errors2 );

                //iterate through all clusters to see if they are in track
                for ( itC = clusterLog_.begin(); itC != clusterLog_.end(); ++itC ) {
                    
                    std::vector< double> point, errors;
                    (itC->second).getPoint( point , errors );
                    
                    double dist = distToLine( point1 , point2 , point );

                    if ( dist < trackRadius_ ) {
                        ctrack_mipids.push_back( itC->first );
                    }

                } //iterate through all clusters to see if they are in track (itC)
                
                //check if plausible track
                if ( ctrack_mipids.size() > minNumClusters_ ) {
                    //create fit for ctrack
                    HcalMipTrack ctrack;
                    for ( std::vector< unsigned int >::iterator it = ctrack_mipids.begin();
                        it != ctrack_mipids.end(); ++it ) {
                        
                        MipCluster* cmip = &clusterLog_[ *it ];
                        for ( int i = 0; i < cmip->getNumHits(); i++ ) {
                            ctrack->addHit( cmip->getHit( iH ) );
                        }//iterate through hits in cluster

                    } //add clusters with mipids to ctrack
                    
                    if ( compMipTracks( best_track , ctrack ) ) {
                        best_track = ctrack;
                        track_mipids = ctrack_mipids;
                        success = true;
                    }//ctrack is better than best_track

                }//ctrack is a plausible track

            } //go through remaining hits as second end point (itC2)
        } //go through all hits as first end point (itC1)
        
        return success;
    }
}

DECLARE_PRODUCER_NS(ldmx, HcalMipTrackProducer);
