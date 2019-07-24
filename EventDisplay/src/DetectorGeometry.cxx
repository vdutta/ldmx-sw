#include "EventDisplay/DetectorGeometry.h"

namespace ldmx {

    DetectorGeometry::DetectorGeometry() {

        hcal_ = new TEveElementList("HCAL");
        sidehcal_ = new TEveElementList("Side HCAL");
        ecal_ = new TEveElementList("ECAL");
        recoilTracker_ = new TEveElementList("Recoil Tracker");
        detector_ = new TEveElementList("LDMX Detector");

        shapeDrawer_ = new EveShapeDrawer();

        drawECAL();
        drawHCAL();
        drawRecoilTracker();

    }

    void DetectorGeometry::drawECAL() {


        for (int col = 0; col < towerXPos.size(); ++col) {
    
            TString colName;
            colName.Form("Tower %d", col);
            TEveGeoShape* hexCol = shapeDrawer_->drawHexPrism(towerXPos[col], towerYPos[col], ecal_front_z+ecal_z_length/2, 0, 0, 0, ecal_z_length, ecal_xy/(3*sqrt(3)), kBlue, 90, colName);

            ecal_->AddElement(hexCol);
        }
    
        detector_->AddElement(ecal_);
    }

    void DetectorGeometry::drawHCAL() {
    
        TEveGeoShape* backHcal = shapeDrawer_->drawRectPrism(
                HCAL_DETECTOR_GEOMETRY.getBoundingBox( HcalSection::BACK ),
                0, 0, 0, kCyan, 90, "Back HCal"); 
        hcal_->AddElement(backHcal);
    
        TEveGeoShape* sideTopHcal = shapeDrawer_->drawRectPrism(
                HCAL_DETECTOR_GEOMETRY.getBoundingBox( HcalSection::TOP ),
                0, 0, 0, kCyan, 90, "Module 1");

        sidehcal_->AddElement(sideTopHcal);
    
        TEveGeoShape* sideBottomHcal = shapeDrawer_->drawRectPrism(
                HCAL_DETECTOR_GEOMETRY.getBoundingBox( HcalSection::BOTTOM ),
                0, 0, 0, kCyan, 90, "Module 4");

        sidehcal_->AddElement(sideBottomHcal);
    
        TEveGeoShape* sideLeftHcal = shapeDrawer_->drawRectPrism(
                HCAL_DETECTOR_GEOMETRY.getBoundingBox( HcalSection::LEFT ),
                0, 0, 0, kCyan, 90, "Module 2");

        sidehcal_->AddElement(sideLeftHcal);
    
        TEveGeoShape* sideRightHcal = shapeDrawer_->drawRectPrism(
                HCAL_DETECTOR_GEOMETRY.getBoundingBox( HcalSection::RIGHT ),
                0, 0, 0, kCyan, 90, "Module 3");

        sidehcal_->AddElement(sideRightHcal);
        hcal_->AddElement(sidehcal_);
    
        detector_->AddElement(hcal_);
    }
    
    void DetectorGeometry::drawRecoilTracker() {
    
        for (int j = 0; j < 4; ++j) {
    
            TString nfront;
            nfront.Form("Stereo%d_front", j+1);
    
            TString nback;
            nback.Form("Stereo%d_back", j+1);
    
            TEveGeoShape *front = shapeDrawer_->drawRectPrism(0, 0, recoilLayerZPos[j]-mono_sep/2, stereo_x_width, stereo_y_width, recoil_sensor_thick, 0, 0, 0, kRed-10, 90, nfront);

    
            if (j % 2 == 0) { // Alternate angle for back layer of a stereo pair.
                TEveGeoShape *back = shapeDrawer_->drawRectPrism(0, 0, recoilLayerZPos[j]+mono_sep/2, stereo_x_width, stereo_y_width, recoil_sensor_thick, 0, 0, stereo_angle, kRed-10, 90, nback);

                recoilTracker_->AddElement(back);
            } else {
                TEveGeoShape *back = shapeDrawer_->drawRectPrism(0, 0, recoilLayerZPos[j]+mono_sep/2, stereo_x_width, stereo_y_width, recoil_sensor_thick, 0, 0, -stereo_angle, kRed-10, 90, nback);

                recoilTracker_->AddElement(back);
            }
    
            recoilTracker_->AddElement(front);
        }
    
        int module1 = 1;
        for (int x = 0; x < 5; ++x) {
            for (int y = 0; y < 2; ++y) {
    
                TString name;
                name.Form("Mono1_%d",module1);
                ++module1;
    
                if (x % 2 != 0) { // Alternate mono layer z by defined separation.
                    TEveGeoShape *front = shapeDrawer_->drawRectPrism(monoSensorXPos[x], monoSensorYPos[y], recoilLayerZPos[4]-mono_sep/2, mono_x_width, mono_y_width, recoil_sensor_thick, 0, 0, 0, kRed-10, 80, name);

                    recoilTracker_->AddElement(front);
                } else {
                    TEveGeoShape *back = shapeDrawer_->drawRectPrism(monoSensorXPos[x], monoSensorYPos[y], recoilLayerZPos[4]+mono_sep/2, mono_x_width, mono_y_width, recoil_sensor_thick, 0, 0, 0, kRed-10, 80, name);

                    recoilTracker_->AddElement(back);
                }
            }
        }
    
        int module2 = 1;
        for (int x = 0; x < 5; ++x) {
            for (int y = 0; y < 2; ++y) {
    
                TString name;
                name.Form("Mono2_%d",module2);
                module2++;
    
                if (x % 2 != 0) { // Alternate mono layer z by defined separation.
                    TEveGeoShape *front = shapeDrawer_->drawRectPrism(monoSensorXPos[x], monoSensorYPos[y], recoilLayerZPos[5]-mono_sep/2, mono_x_width, mono_y_width, recoil_sensor_thick, 0, 0, 0, kRed-10, 90, name);

                    recoilTracker_->AddElement(front);
                } else {
                   TEveGeoShape *back = shapeDrawer_->drawRectPrism(monoSensorXPos[x], monoSensorYPos[y], recoilLayerZPos[5]+mono_sep/2, mono_x_width, mono_y_width, recoil_sensor_thick, 0, 0, 0, kRed-10, 90, name);
                    recoilTracker_->AddElement(back);
                }
            }
        }
    
        detector_->AddElement(recoilTracker_);
    }
}
