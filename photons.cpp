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
    //vector <float> event_position;

    // experimenting with using ROOT 3D physics vectors instead of C++
    TVector3 event_position; 

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
            event_position.SetXYZ(R*sin(theta)*cos(phi), R*sin(theta)*sin(phi), R*cos(theta));

            // print the vector to the terminal nicely for sanity checking 
            printf("The vector [x, y, z] is: [");
            
            // empty array with 3 elements 
            double d[3];

            // fill d-array with coordinates of event
            //event_position.GetCoordinates(d);
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
            //ROOT::Math::Polar3DVector intercepts[num_photons]; 
            float time_of_flight[num_photons]; 
            gammaIntercept(time_of_flight); 

            puts("INTERCEPTS AT: \n");
            for (int i =0; i < num_photons; i++){
                float temp[3];
                //intercepts[i].GetCoordinates(temp);
                printf("intercept %d at: [ %f, %f, %f ]\n", i, intercepts[i][0], intercepts[i][1], intercepts[i][2]);
            }

            // checking 
            // printf("returning the origin within the event\n");
            // double *origin = get_origin(); 
            // for (int i = 0; i < 3; i++){
            //     printf("\nloopin #%d\n", i);
            //     printf("%f, ", *origin);
            //     ++origin;

            // };
            // printf("\nFINISH\n"); 
        };

        int convertGammas(){
            // convert to a number of photons given the event energy
            int photons_per_mev = 5;  

            const int num_photons = energy * photons_per_mev; 
            //const int num_photons = 10;
            printf("Number of photons is: %d\n", num_photons);

            return num_photons;
        };

        void gammaDir(){
            //vector <float> gamma_direction; 
            TVector3 gamma_direction; 
            TRandom *r1 = new TRandom(time(NULL));

            for(int i =0; i < num_photons; i++){
                // create random r, theta, phi relative to event location
                // float r     = 1.0; // doesn't really matter for the direction vector
                // float theta = rand() / (static_cast <float> (RAND_MAX) / M_PI);
                // float phi   = rand() / (static_cast <float> (RAND_MAX) / (2 * M_PI)); 
                
                float phi = r1->Uniform(0, 2*pi);
                float theta = r1->Uniform(0,pi);
                float r = r1->Uniform(0,6);
                // populate vector and normalise
                // float vec_size = sqrt (pow(r,2) + pow(theta,2) + pow(phi,2));
                // gamma_direction = {r/vec_size, theta/vec_size, phi/vec_size};

                gamma_direction.SetXYZ(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta));
                gamma_direction = gamma_direction.Unit(); 

                // add vector to the containing array for the event 
                gamma_dirs.push_back(gamma_direction);
            };
        }; 

        void gammaIntercept(float *time_of_flights){
            // calculates the intercept of photon and PMTs
            /* Takes a POINTER to the empty fixed-size C-array container interceptsand time of flights. By default, the pointer 
            points to the first element of intercepts. The pointer itself is incremented at the end of the 
            for loop, so it's pointing at the next element of intercepts for assignment. The maths for the intercept is 
            taken from wikipedia sphere/line intercept page. */

            TVector3 intercept;  // intercept of 1 photon with PMT sphere surface 
            TVector3 centre;     // centre of sphere (detector) and coordinate system
            TVector3 diff;       // vector difference between event origin and sphere centre 
            float d;                              // distance along line from event origin to intercept  
            float tof;                            // time of flight 
            float del;                            // huge chunk of equation I'd previously forgotten to add...

            centre.SetXYZ(0,0,0);
            // loop over all the photons 
            int flag = 1; 
            for (auto i : gamma_dirs){
                diff = event_position - centre;
                
                del = pow((i.Dot(diff)), 2) - (diff.Mag2() - 36);
                if (flag == 1){
                    cout << "dir vec is: "<< endl;
                    for (int k = 0; k < 3; k++){
                        cout << "[";
                        cout << i[k] << ",";

                    };
                    cout << endl; 
                    cout << "Magnitude of allegeded unit vector is: " << i.Mag() << endl;
                    flag = 2; 
                };
                float d1 = - (i.Dot(diff)) + sqrt(del);
                float d2 = - (i.Dot(diff)) - sqrt(del); 
                printf("del is: %f\nd1: %f \nd2: %f", del, d1, d2);

                double checker = (event_position + d1*i -centre).Mag();
                printf("\nChecker is %f", checker); 
                if (checker == 6){
                    d = d1;
                    puts("\n USING d1!\n");
                } else{
                    d = d2; 
                    puts("\n USING d2!\n");

                };
                 
                // intercept is event_origin + ( unit dir vec * d )
                TVector3 intercept1 = event_position + (d1 * i);  
                TVector3 intercept2 = event_position + (d2 * i); 

                cout << "\nINT1 r = " << sqrt(pow(intercept1[0],2) + pow(intercept1[1],2) +pow(intercept1[2], 2)) << endl; 
                cout << "INT2 r = " << sqrt(pow(intercept2[0],2) + pow(intercept2[1],2) +pow(intercept2[2], 2)) << endl; 
                //printf("\nint1 is: %f\nint2 is: %f ", intercept1[0], intercept2[0]);
                intercepts.push_back(intercept1); 
                //intercepts.push_back(intercept2); 
                //intercept[0] = 6;
                

                //printf("Distance to sphere surface: %f\n", d); 

                // calculate time of flight and assign to correct array element
                tof = d / pow(3,8); 
                *time_of_flights = tof; 

                
                

                // print out the intercept 
                // float d[3];
                // intercept.GetCoordinates(d);
                // cout << "Intercept: "; 
                // for (auto j : d){
                //     cout << j << " ";
                // };
                // cout << endl;  

                // add intercept of photon to intercepts list via pointer
                //*intercepts = intercept; 

                // increment pointer to point at next element of intercepts
                //++intercepts; 
                ++time_of_flights;   
            };
        };

        double * get_origin(){

            // MUST BE STATIC for this to work!
            static double x[3];
            //event_position.GetCoordinates(x); 

            // puts("INSIDE GET ORIGIN!");
            for (int i = 0; i < 3; i++){
                x[i] = event_position[i];
            }; 
            
            return x;    
        };

        double * get_intercepts(int number){

            // weird...pass the intercept number in a loop in sphere, and this function 
            // will return a pointer to an array containing the intercept of that photon
            //puts("INTERCEPT AT: \n");
            static double temp[3];
            //intercepts[number].GetCoordinates(temp);
            //printf("intercept %d at: [ %f, %f, %f ]\n", i, temp[0], temp[1], temp[2]);
            for (int i = 0; i < 3; i++){
                temp[i] = intercepts[number][i];
            };
            return temp; 
        };
        
        int get_num_photons(){
            return num_photons;
        };
    };