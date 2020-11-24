
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVirtualGeoTrack.h"
#include "./photons.cpp"
#include <math.h>
#include <array>
using namespace std; 

void plot_event(int num_events){
    
    // initialise a list of Event objects of a given number
    Event *events = new Event[num_events];
    for (int i = 0; i < num_events; i++){
        double * origin = events[i].get_origin();
        int num_photons =  events[i].get_num_photons(); 
        cout << "NUMBER OF INTERCEPT POINTS RETURNED " << num_photons << endl; 

        // container for all the photon-sphere intercepts
        double cart1[num_photons][3]; 

        // print the photon-sphere intercepts 
        cout << "The photon intercepts are: " << endl;

        // loop over each photon-sphere intercept
        for (int l = 0; l < num_photons; l++){
        
            // obtain singlular photon-sphere intercept
            double * intercept = events[i].get_intercepts(l); 
        
            // populate container with intercepts
            cart1[l][0] = intercept[0]; 
            cart1[l][1] = intercept[1]; 
            cart1[l][2] = intercept[2];

            // print the photon-sphere intercepts 
            cout << "[";
            for (int k =0; k < 3; k++){
                cout << cart1[l][k];
                if (k == 2){
                    cout << "]" << endl;
                } else {
                    cout << ", ";
                };
            }; 
        };

        cout << "Event Origin passed: ";
        for (int j = 0; j < 3; j++){
            cout << origin[j] << " ";
        };
        cout << endl;

        // loop over every photon intercept
        // create a pointer to an empty array that will store all of the TPoly3D lines 
        TPolyLine3D *trajectories[num_photons];
        for (int d = 0; d < num_photons; d++){
            // array stores the coordinates of 1 point
            double coords[3]; 

            // add a line to the array of lones 
            trajectories[d] = new TPolyLine3D(2);
            TPolyMarker3D *point = new TPolyMarker3D(1, cart1[d], 3);
            trajectories[d]->SetPoint(0, origin[0], origin[1], origin[2]);
            trajectories[d]->SetPoint(1, cart1[d][0], cart1[d][1], cart1[d][2]); 
            trajectories[d]->SetLineColor(i+1);
            point->SetMarkerColor(i+1);
            // draw the line 
            trajectories[d]->Draw();
            point->Draw();
        };

        // create a PolyPoint to mark event origin 
        TPolyMarker3D *pm3d2 = new TPolyMarker3D(1, origin, 5);
        pm3d2->Draw();
        };

    // garbage collection 
    delete[] events;
};

void sphere(){

    // create a canvas 
    TCanvas *c1 = new TCanvas("c1", "Electron Event Visualisation", 200, 10, 700, 500);
    
    // create pad within the canvas 
    TPad *p1 = new TPad("p1", "p1", 0.05, 0.02, 0.95, 0.82, 46, 3, 1); 
    
    // geometry manager for sphere
    gSystem->Load("libGeom");
    new TGeoManager("world", "the simplest geometry");

    // these don't matter for this purpose but need to be defined 
    TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
    TGeoMedium   *med = new TGeoMedium("Vacuum",1,mat);

    // create sphere of defined medium with r[0->6], theta[0->180], phi[0->360]
    TGeoVolume *top=gGeoManager->MakeSphere("Top",med,0.,6,0, 180., 0., 360.);
    gGeoManager->SetTopVolume(top);

    // always close geomtery after creation to allow it to be drawn
    gGeoManager->CloseGeometry();
    top->SetLineColor(kMagenta);

    // default to invisible, so set sphere visible  
    gGeoManager->SetTopVisible();

    // draw sphere
    top->Draw();

    // now GENERATE EVENTS 
    int num_events = 3;
    plot_event(num_events); 
};

int main(){ 
    sphere();
}; 