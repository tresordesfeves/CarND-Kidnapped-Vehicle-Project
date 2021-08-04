/* particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 ./clean.sh
 ./build.sh.
 ./run.sh
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::max;



  std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */

// This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(x, std[0]);
// This line creates a normal (Gaussian) distribution for y
  normal_distribution<double> dist_y(y, std[1]);
// This line creates a normal (Gaussian) distribution for theta
  normal_distribution<double> dist_theta(theta, std[2]);

// TODO: Set the number of particles
  num_particles = 100; 
  heaviest_particle=999;
// create a shell variable "particle" using the type " struct Particle" 
  Particle particle;


  for (int i = 0; i < num_particles; ++i) {

    
    particle.id=i; // tag this very particle
    
    //  Sample x from a normal distribution using the random engine gen initialized earlier.
    particle.x= dist_x(gen);
    //  Sample y from a normal distribution using the random engine gen initialized earlier.
    particle.y= dist_y(gen);
    //  Sample theta from a normal distribution using the random engine gen initialized earlier.
    particle.theta= dist_theta(gen);
    // Initialize weight to 1
    particle.weight=1;

    // Update the list of particles (particle vector)
    particles.push_back(particle);

   
  }




}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */


  //!!!! double sigma_pos [3] = {0.3, 0.3, 0.01};



    /*  we add gaussian noise for each dimension of the particles positions
      Standard deviation is set for each position dimension (x,y,theta ) once and for all particles
    */

    normal_distribution<double> dist_x(0, std_pos[0]); //define the gaussian distributin for x
    normal_distribution<double> dist_y(0, std_pos[1]); //define the gaussian distributin for y
    normal_distribution<double> dist_theta(0, std_pos[2]); //define the gaussian distributin for theta

    double drift,delta_theta,new_theta;

if (abs(yaw_rate) < 1.0e-5) // only test once per time step instead of per particle
  {
    for (int i = 0; i < num_particles; ++i)
      {
        particles[i].x+= velocity * cos(particles[i].theta)* delta_t;
        particles[i].y+= velocity * sin(particles[i].theta)* delta_t;

      // generate and add the ramdom noise to each [article position:
        particles[i].x+=dist_x(gen);
        particles[i].y+=dist_y(gen);
        particles[i].theta+=dist_theta(gen);

      }
  }
else
  {
      drift = velocity/yaw_rate;
      delta_theta=yaw_rate* delta_t;// calculate the variation of direction for all particles for this time step 
  
      for (int i = 0; i < num_particles; ++i)
      //*********************
      // generate the noise for x,y,theta ( normal dist ( mu =x,y,theta), sigma= std )
      // or should only be for each update like here 
        {      
            // update each particle position withthe motion:   
              new_theta+= delta_theta+ particles[i].theta ; //update particule direction 
              particles[i].x += drift*(sin(new_theta)-sin(particles[i].theta));
              particles[i].y += drift*(cos(particles[i].theta)-cos(new_theta));
              particles[i].theta = new_theta;

            // generate and add the ramdom noise to each [article position:
              particles[i].x+=dist_x(gen);
              particles[i].y+=dist_y(gen);
              particles[i].theta+=dist_theta(gen);
        }
    }



  }

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  int predicted_count = predicted.size();
  int observations_count = observations.size();
  int closest_LM_index;

  double distance_obs_LM;

  // this is for a specific particle

  for ( int i = 0; i < observations_count ; i++) // for each radar observation (ie : observed measurement) 
    { 
      double shortest_distance= 999999; // initialize for each measurement


      for ( int j = 0; j < predicted_count; ++j) // for each map landmark within the range of this specific particle 
        {
         distance_obs_LM= dist(predicted[j].x,predicted[j].y, observations[i].x,observations[i].y);

         if (distance_obs_LM<shortest_distance)
            {
              shortest_distance=distance_obs_LM;
              closest_LM_index=j ;//track the index within the in-range landmarks vector    

            }
        } 

        observations[i].id=closest_LM_index;
    } 

  }

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                    const vector<LandmarkObs> &observations, 
                                    const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
    int map_size=map_landmarks.landmark_list.size(); 
    int obs_size = observations.size();
    double xp,yp,xc,yc,xm,ym; 
    double theta, cos_theta,sin_theta;
    double mu_x,  mu_y;
    double weight;
    double total_weight=0;
    LandmarkObs LM; // placeholder for calculations 
     
    vector<LandmarkObs> LM_in_particle_i_range; // landmarks within the particle sensor range
    vector<LandmarkObs> particle_i_observations; // obvervations once applied to a specific particle and converted in map coordinates


for (int i = 0; i < num_particles; ++i) // FOR EACH PARTICLE  !!!!!!
  { /* 
    in regard to this very particle i , predict measurements to all the map landmarks within the sensor range 
     and transdorm the prediction in map coordinates: 
    */

    weight=1;

    xp = particles[i].x ; 
    yp = particles[i].y ; 
    theta = particles[i].theta ; 
    cos_theta =cos(theta),sin_theta = sin(theta);// initialized ar the particle level to avoid recalculation at each radar observation
    // for this particle i : predict measurements to all the map landmarks within sensor range 

    for (int j = 0; j < map_size; ++j) // predict all landmarks within the range of particle i 



        {
         
          if (dist( xp,  yp,  map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f) < sensor_range )
          { 

            /* although both sructures are semanticaly identical, the field names are different 
            we stick to one type for conveniency in :
            ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) */ 

            LM.id=map_landmarks.landmark_list[j].id_i;
            LM.x=map_landmarks.landmark_list[j].x_f;
            LM.y=map_landmarks.landmark_list[j].y_f;

            LM_in_particle_i_range.push_back(LM);
          }
            //LM_in_particle_i_range.push_back(map_landmarks.landmark_list[j])}

        }

    // all measurements from the radar in RESPECT TO THIS VERY PARTICLE converted in map coordinate: 


    for (int j = 0; j < obs_size; ++j) // for each observation
      {
      //(xc,yc) car radar observations in the car coordinate system
      //(xp,yp) particle cordinate in the map system
      //(xm,ym) observations measurements applied to the particle and transformed to the map system

        xc=observations[j].x;
        yc=observations[j].y;

        xm= xp + (cos_theta * xc ) - ( sin_theta * yc);
        ym= yp + (sin_theta * xc ) + ( cos_theta * yc); 

        LM.x=xm;
        LM.y=ym;
        particle_i_observations.push_back(LM);

        }

      dataAssociation(LM_in_particle_i_range,particle_i_observations);


/*
      calculate the weight for this very particle

      multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y)
*/

      // all info is in particle_i_observations: now calculate the weight


      for (int j = 0; j < obs_size; ++j)
        {
          mu_x = LM_in_particle_i_range[particle_i_observations[j].id].x;
          mu_y = LM_in_particle_i_range[particle_i_observations[j].id].y;
          
          weight *=multiv_prob(std_landmark[0],std_landmark[1], particle_i_observations[j].x , particle_i_observations[j].y, mu_x, mu_y);

        }

      particles[i].weight=weight; // update the weight of particle i
      total_weight+=weight; // track total weight of all particles 
      LM_in_particle_i_range.clear();// start fresh on the next particle
      particle_i_observations.clear();// start fresh on the next particle

    }

  heaviest_particle=0;

  for (int i = 0; i < num_particles; ++i)
    {
      particles[i].weight=particles[i].weight/total_weight;// normalize to probability distribution across all particles

      heaviest_particle= max(particles[i].weight , heaviest_particle); //update the highest weight ( needed  at the next resampling step )

    }

  }



void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  std::vector<Particle> resampled_particules; // to collect the new set of particles 

  std::uniform_int_distribution<int> int_distribution(0,num_particles-1); 


  double beta =0;
  std:: uniform_real_distribution<double> real_distribution(0.0, 2* heaviest_particle);

  int index = int_distribution(gen); // random integer within the O to particle count range  


 for (int i=0; i<num_particles-1; ++i) {
  
    beta += real_distribution(gen) ; // random real within the range of O to twice the heaviest particle range; 
    while(beta>particles[index].weight) 

      {
        beta -= particles[index].weight;// lighter particles have smaller circonference portion , less change to be picked ( and the opposite for heavier one  )  
  
        index=(index+1)%num_particles;

      }
    
    resampled_particules.push_back(particles[index]);
    }

  particles=resampled_particules;
   /* sending out the same number of particles:heavier ones have more chance to be duplicated 
   lighters ones to be dropped */


}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

