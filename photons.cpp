#include <string>
#include <cstdlib> 
#include <vector>
#include <ctime>
#include <math.h>
#include <cmath>
#include <iostream>
#include "TCanvas.h"
#include "TVector3.h"
#include "Math/Vector3D.h"
#include "TRandom.h"
#include "TMath.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TStyle.h"
#include "TSystem.h"
using namespace std; 


class Event {
    // private data members
    // Events have an energy and a position vector 
    float energy;
    int num_photons; 
    float pi = TMath::Pi();
    vector < ROOT::Math::Polar3DVector > gamma_dirs; // array of vectors, each row is a polar 3D vector  
    
    //vector <float> event_position;

    // experimenting with using ROOT 3D physics vectors instead of C++
    ROOT::Math::Polar3DVector event_position; 

    // constructor of class randomly populates the energy and position attributes, converts E --> # gammas and 
    // specifies a random direction for the photons 
    public: 
        Event(){
            // initialise the random number generator using the system time seed so it's different each time 
            //srand(time(NULL));

            // randomly assign an energy between 1 --> 5MeV 
            // rand() is defined in <cstdlib> and generates a random number from 0 --> RAND_MAX constant 
            // taking RAND_MAX % 5 is number from 0 --> 4, so add 1 for 1 --> 5. 
            //energy = 1 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX) / 4);   

            // ROOT random number generator 
            // define a new TRandom class type via a pointer r1 
            TRandom *r1 = new TRandom(time(NULL));

            // get energy from uniform distribution between 0-->5MeV
            energy = r1->Uniform(1, 5); 
            
            printf("The energy from ROOT is: %f\n", energy);
            // generate random phi between 0 --> 2pi
            // generate random theta between 0 --> pi  
            // generate random R between 0 --> 6
            // float phi   = rand() / (static_cast <float> (RAND_MAX) / (2 * 3.141592653589793238463)); 
            // float theta = rand() / (static_cast <float> (RAND_MAX) / 3.141592653589793238463); 
            // float R     = rand() / (static_cast <float> (RAND_MAX) / 6); 
            float phi = r1->Uniform(0, 2*pi);
            float theta = r1->Uniform(0,pi);
            float R = r1->Uniform(0,6);
            
            // create the spherical coordinate position vector (R, theta, phi)
            event_position.SetCoordinates(R, theta, phi);

            // print the vector to the terminal nicely for sanity checking 
            printf("The vector [R, theta, phi] is: [");
            
            // empty array with 3 elements 
            double d[3];

            // fill d-array with coordinates of event
            event_position.GetCoordinates(d);
            for (int i = 0; i < 3; i++){
                if (i < 2){
                    printf("%f ", d[i]);
                } else{
                    printf("%f", d[i]);
                }
            };
            printf("]\n"); 
            
            // converts event_energy to the number of photons 
            convertGammas(); 

            // assign a direction to each gamma
            gammaDir();

            // print all the unit photon direction vectors 
            // for (auto i : gamma_dirs){
            //     double d[3]; 
            //     i.GetCoordinates(d);
            //     for (auto j : d){
            //         cout << j << " ";
            //     };
            //     cout << endl; 
            // };

            // find intercept of photon direction and PMT sphere 
            // create an array of ROOT 3D polar vectors which store the photon/PMT intercepts
            // and calculate the time of flight and store in a C-array 
            ROOT::Math::Polar3DVector intercepts[num_photons]; 
            float time_of_flight[num_photons]; 
            gammaIntercept(intercepts, time_of_flight); 
        };

        void convertGammas(){
            // convert to a number of photons given the event energy
            int photons_per_mev = 100;  

            num_photons = energy * photons_per_mev; 
            printf("Number of photons is: %d\n", num_photons);
        };

        void gammaDir(){
            //vector <float> gamma_direction; 
            ROOT::Math::Polar3DVector gamma_direction; 
            TRandom *r1 = new TRandom(time(NULL));

            for(int i =0; i < num_photons; i++){
                // create random r, theta, phi relative to event location
                // float r     = 1.0; // doesn't really matter for the direction vector
                // float theta = rand() / (static_cast <float> (RAND_MAX) / M_PI);
                // float phi   = rand() / (static_cast <float> (RAND_MAX) / (2 * M_PI)); 
                
                float phi = r1->Uniform(0, 2*pi);
                float theta = r1->Uniform(0,pi);
                
                // populate vector and normalise
                // float vec_size = sqrt (pow(r,2) + pow(theta,2) + pow(phi,2));
                // gamma_direction = {r/vec_size, theta/vec_size, phi/vec_size};

                gamma_direction.SetCoordinates(1.0, theta, phi);
                gamma_direction = gamma_direction.unit(); 

                // add vector to the containing array for the event 
                gamma_dirs.push_back(gamma_direction);
            };
        }; 
        void gammaIntercept(ROOT::Math::Polar3DVector *intercepts, float *time_of_flights){
            // calculates the intercept of photon and PMTs
            /* Takes a POINTER to the empty fixed-size C-array container interceptsand time of flights. By default, the pointer 
            points to the first element of intercepts. The pointer itself is incremented at the end of the 
            for loop, so it's pointing at the next element of intercepts for assignment. The maths for the intercept is 
            taken from wikipedia sphere/line intercept page. */

            ROOT::Math::Polar3DVector intercept;  // intercept of 1 photon with PMT sphere surface 
            ROOT::Math::Polar3DVector centre;     // centre of sphere (detector) and coordinate system
            ROOT::Math::Polar3DVector diff;       // vector difference between event origin and sphere centre 
            float d;                              // distance along line from event origin to intercept  
            float tof;                            // time of flight 
            
            // loop over all the photons 
            for (auto i : gamma_dirs){
                diff = event_position - centre;
                d = - (i.Dot(diff));

                // calculate time of flight and assign to correct array element
                tof = d / pow(3,8); 
                *time_of_flights = tof; 

                // intercept is event_origin + ( unit dir vec * d )
                intercept = event_position + (d * i);

                // print out the intercept 
                // float d[3];
                // intercept.GetCoordinates(d);
                // cout << "Intercept: "; 
                // for (auto j : d){
                //     cout << j << " ";
                // };
                // cout << endl;  

                // add intercept of photon to intercepts list via pointer
                *intercepts = intercept; 

                // increment pointer to point at next element of intercepts
                ++intercepts; 
                ++time_of_flights;   
            };
        };
        
        void visualise_event(){
            /* Creates a 3D graphical representation of the event location, photon paths and intercepts using ROOTs
            3D visualisations. */
            
            // load geometry library 
            gSystem->Load("libGeom");

            // create instance of geometry manager class
            new TGeoManager("world", "the simplest geometry");

            // volume material and medium - can be ignored for now 
            TGeoMaterial *mat = new TGeoMaterial("Vacuum", 0,0,0);
            TGeoMedium *med = new TGeoMedium("Vacuum", 0,0,0);  

            // define rmin and rmax of sphere in units of cm
            Double_t rmin = 0; 
            Double_t rmax = 600;

            // make sphere, args are [name, medium, rmin, rmax, theta min, theta nax, phi min, phi max]
            TGeoVolume *top = gGeoManager->MakeSphere("Top", med, rmin, rmax, 0, 180, 0, 360);

            // make this volume our "world"
            gGeoManager->SetTopVolume(top);
            gGeoManager->CloseGeometry(); // some random optimisation checks 

            // display the world 
            top->SetLineColor(kMagenta);
            gGeoManager->SetTopVisible(); // default top geometry invisible
            top->Draw(); 
        };
    }; 


int main () {
    puts("Hello world\n");

    // create an event object 
    Event ev1; 

    // access visualisation
    ev1.visualise_event();

}
