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
#include "TRandom3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TApplication.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVector3.h"
using namespace std; 

class Event {
    // private data members
    // Events have an energy and a position vector 
    float energy;
    int num_photons; 
    float pi = TMath::Pi();
    vector < TVector3 > gamma_dirs; // array of vectors, each element is a polar 3D vector  
    vector < TVector3 > intercepts; // array of vectors, each element is photon/sphere intercept coords
    vector <float> time_of_flight;  // vector stores time of flight information 

    // experimenting with using ROOT 3D physics vectors instead of C++
    TVector3 event_position; 

    // constructor of class randomly populates the energy and position attributes, converts E --> # gammas and 
    // specifies a random direction for the photons 
    public: 
        Event(){ 
            // ROOT random number generator 
            // define a new TRandom class type via a pointer r1
            // zero seed - random number each time called  
            TRandom3 *r1 = new TRandom3(0);

            // get energy from uniform distribution between 0-->5MeV
            energy = r1->Uniform(1, 5); 
            
            printf("The energy from ROOT is: %f\n", energy);
            
            // randomly produce the event position coordinates in spherical polars
            float phi = r1->Uniform(0, 2*pi);
            float theta = r1->Uniform(0,pi);
            float R = r1->Uniform(0,6);
            
            // create the spherical coordinate position vector (R, theta, phi)
            event_position.SetXYZ(R*sin(theta)*cos(phi), R*sin(theta)*sin(phi), R*cos(theta));

            // print the vector to the terminal nicely for sanity checking 
            printf("The vector [x, y, z] is: [");
            for (int i = 0; i < 3; i++){
                if (i < 2){
                    printf("%f ", event_position[i]);
                } else{
                    printf("%f", event_position[i]);
                }
            };
            printf("]\n"); 
            
            // converts event_energy to the number of photons 
            num_photons = convertGammas(); 

            // assign a direction to each gamma
            gammaDir();

            // find intercept of photon direction and PMT sphere 
            // and calculate the time of flight and store in a vector 
            gammaIntercept(); 

            // print out the calculated intercepts 
            puts("INTERCEPTS AT: \n");
            for (int i =0; i < num_photons; i++){
                printf("intercept %d at: [ %f, %f, %f ]\n", i, intercepts[i][0], intercepts[i][1], intercepts[i][2]);
            };
        };

        int convertGammas(){
            // convert to a number of photons given the event energy
            int photons_per_mev = 5;  
            const int num_photons = energy * photons_per_mev; 
            printf("Number of photons is: %d\n", num_photons);
            
            return num_photons;
        };

        void gammaDir(){
            TVector3 gamma_direction; 
            TRandom *r1 = new TRandom(0);

            // for each photon generate a direction vector
            for(int i =0; i < num_photons; i++){
                float phi = r1->Uniform(0, 2*pi);
                float theta = r1->Uniform(0,pi);
                float r = r1->Uniform(0,6);

                // convert from spherical --> cartesian coords for easy time later
                gamma_direction.SetXYZ(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta));

                // convert to unit vector 
                gamma_direction = gamma_direction.Unit(); 

                // add vector to the containing vector for the event 
                gamma_dirs.push_back(gamma_direction);
            };
        }; 

        void gammaIntercept(){
            // calculates the intercept of photon and PMTs
            /* Takes a POINTER to the empty fixed-size C-array container interceptsand time of flights. By default, the pointer 
            points to the first element of intercepts. The pointer itself is incremented at the end of the 
            for loop, so it's pointing at the next element of intercepts for assignment. The maths for the intercept is 
            taken from wikipedia sphere/line intercept page. */

            TVector3 intercept;   // intercept of 1 photon with PMT sphere surface 
            TVector3 centre;      // centre of sphere (detector) and coordinate system
            TVector3 diff;        // vector difference between event origin and sphere centre 
            float d;              // distance along line from event origin to intercept  
            float tof;            // time of flight 
            float del;            // huge chunk of equation I'd previously forgotten to add...

            centre.SetXYZ(0,0,0); // center of containing sphere is the origin

            // loop over all the photons 
            for (auto i : gamma_dirs){

                // mathematics of intercept from wikipedia 
                diff = event_position - centre;
                del = pow((i.Dot(diff)), 2) - (diff.Mag2() - 36);

                // there are two solutions, a forwards and a backwards one +-del 
                float d1 = - (i.Dot(diff)) + sqrt(del);
                float d2 = - (i.Dot(diff)) - sqrt(del); 

                // take the positive value solution as the correct one
                if (d1 > 0){
                    d = d1;
                } else{
                    d = d2;
                };
                 
                // intercept is event_origin + ( unit dir vec * d )
                TVector3 intercept = event_position + (d * i);  
                intercepts.push_back(intercept); 
                
                // calculate time of flight and assign to correct array element
                tof = d / pow(3,8); 
                time_of_flight.push_back(tof); 
            };
        };

        double * get_origin(){
            // MUST BE STATIC for this to work!
            static double x[3];

            // puts("INSIDE GET ORIGIN!");
            for (int i = 0; i < 3; i++){
                x[i] = event_position[i];
            }; 
            
            return x;    
        };

        double * get_intercepts(int number){
            // weird...pass the intercept number in a loop in sphere, and this function 
            // will return a pointer to an array containing the intercept of that photon
            static double temp[3];
            
            for (int i = 0; i < 3; i++){
                temp[i] = intercepts[number][i];
            };
            return temp; 
        };
        
        int get_num_photons(){
            return num_photons;
        };
    };