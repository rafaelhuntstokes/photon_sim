
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

void sphere(){

    // create a canvas 
    TCanvas *c1 = new TCanvas("c1", "Electron Event Visualisation", 200, 10, 700, 500);
    
    // create pad within the canvas 
    TPad *p1 = new TPad("p1", "p1", 0.05, 0.02, 0.95, 0.82, 46, 3, 1); 
    
    // geometry manager for sphere
    gSystem->Load("libGeom");
    new TGeoManager("world", "the simplest geometry");
    TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
    TGeoMedium   *med = new TGeoMedium("Vacuum",1,mat);
    TGeoVolume *top=gGeoManager->MakeSphere("Top",med,0.,6,0, 180., 0., 360.);
    gGeoManager->SetTopVolume(top);
    gGeoManager->CloseGeometry();
    top->SetLineColor(kMagenta);
    gGeoManager->SetTopVisible();
    top->Draw();

    Event ev1;
    double * origin = ev1.get_origin();
    int num_photons =  ev1.get_num_photons(); 
    cout << "NUMBER OF INTERCEPT POINTS RETURNED " << num_photons << endl; 

    // container for all the photons
    double cart1[num_photons][3]; 
    for (int i = 0; i < num_photons; i++){
        // obtain photon-sphere intercept
        double * intercept = ev1.get_intercepts(i); 
        
        // convert to cartesian 
        cart1[i][0] = intercept[0]; //* sin(intercept[1]) * cos(intercept[2]);
        cart1[i][1] = intercept[1]; //* sin(intercept[1]) * sin(intercept[2]);
        cart1[i][2] = intercept[2]; //* cos(intercept[1]); 

        
    };
    //int s_array = cart1.size();
    //cout << "Number of points plotted " << s_array << endl; 
    cout << "The photon intercepts are: " << endl;
    for (int j =0; j < num_photons;  j++){
        cout << "[";
        for (int k =0; k < 3; k++){
            cout << cart1[j][k] << " "; 
        }; 
        cout << "]"; 
        cout << endl; 
    };

    // convert to cartesian coordinates 
    float cart[3]; 
    cart[0] = origin[0]; //* sin(origin[1]) * cos(origin[2]);
    cart[1] = origin[1]; //* sin(origin[1]) * sin(origin[2]);
    cart[2] = origin[2]; //* cos(origin[1]); 

    cout << "Event Origin passed: ";
    for (int j = 0; j < 3; j++){
        cout << cart[j] << " ";
    };
    cout << endl;

    // loop over every photon intercept

    // create a pointer to an empty array that will store all of the TPoly3D lines 
    TPolyLine3D *trajectories[num_photons];
    for (int i = 0; i < num_photons; i++){
        // array stores the coordinates of 1 point
        double coords[3]; 
        
        // loop over each coordinate
        for (int j = 0; j < 3; j++){
            coords[j] = cart1[i][j]; 
        }; 

        // draw a line from the origin to the intercept
        // float r = sqrt(pow(coords[0], 2) + pow(coords[1], 2) + pow(coords[2], 2));
        // cout << "r is  " << r << endl; 

        // add a line to the array of lones 
        trajectories[i] = new TPolyLine3D(2);
        TPolyMarker3D *point = new TPolyMarker3D(1, coords, 3);
        trajectories[i]->SetPoint(0, cart[0], cart[1], cart[2]);
        trajectories[i]->SetPoint(1, cart1[i][0], cart1[i][1], cart1[i][2]); 
        // pl3d->SetLineWidth(1);
        // pl3d->SetLineColor(2);
        // if (i == 0){
        //     pl3d->SetPoint(1, 1, 2, 3);
        //     pl3d->SetPoint(2, cart[0], cart[1], cart[2]);
        //     cout << cart[0] << " " << cart[1] << " " << cart[2] << endl; 
        // pl3d->Draw();
        // point->Draw();
        // };

        // draw the line 
        trajectories[i]->Draw();
        point->Draw();

        //delete pl3d;
    };



    // // create a PolyPoint to mark event origin 
    TPolyMarker3D *pm3d2 = new TPolyMarker3D(1, cart, 5);
    pm3d2->Draw();
    // // create  a 3D polyline with  
    // TPolyLine3D *pl3d1 = new TPolyLine3D(2, cart);
    // TPolyMarker3D *posPhoton = new TPolyMarker3D(1, cart1, 3);
    

    // // line width and colour attributes 
    // pl3d1->SetLineWidth(1);
    // pl3d1->SetLineColor(2);

    // // draw
    
    // pl3d1->Draw(); 
    
    // posPhoton->Draw(); 

};

int main(){
     
    sphere();
}; 