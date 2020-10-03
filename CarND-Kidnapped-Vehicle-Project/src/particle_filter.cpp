/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
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
using std::normal_distribution;

using std::string;
using std::vector;

#define PARTICLES_NUMBERS 60
#define EPS 0.000001

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
   std::default_random_engine gen;
   num_particles = PARTICLES_NUMBERS;  // TODO: Set the number of particles

    
   normal_distribution<double> dist_x(x, std[0]);
   normal_distribution<double> dist_y(y, std[1]);
   normal_distribution<double> dist_theta(theta, std[2]);

   particles.resize(num_particles); 
   weights.resize(num_particles);
   double init_weight = 1.0;
   for (int i = 0; i < num_particles; i++){
        particles[i].id = i;
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
        particles[i].weight = init_weight;
   }
   is_initialized = true;
     
    // Print your samples to the terminal.
   //std::cout << "Sample " << " " << x << " " << y << " " << theta << std::endl;
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
  
  std::default_random_engine gen;
  num_particles = PARTICLES_NUMBERS;  // TODO: Set the number of particles

    
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  
  double yaw_delta = yaw_rate*delta_t;
  double vel_delta = velocity*delta_t;
  
  for (int i=0;i<num_particles;i++)
  {
    if (yaw_rate < EPS)
    {
      particles[i].x += vel_delta * (cos(particles[i].theta));
      particles[i].y += vel_delta * (sin(particles[i].theta));
    }
    else
    {
      double vel_yaw_rate = velocity/yaw_rate;
      particles[i].x += vel_yaw_rate * (sin(particles[i].theta + yaw_delta) - sin(particles[i].theta));
      particles[i].y += vel_yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_delta));
      particles[i].theta += yaw_delta;
    }
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);   
    particles[i].theta += dist_theta(gen);
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
  
   ////// logic implemented as part of updateWeights function itself

}
                                     
/*double cal_distance(double x1, double y1, double x2, double y2) {
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}*/

                                     
                                     
double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs, double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}
                                     

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your   particles are located according to the MAP'S coordinate system. 
     You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  //weights.clear(); 
  
  for (int i=0;i<num_particles;i++)
  {
    
    double xp = particles[i].x;
    double yp = particles[i].y;
    double thetap = particles[i].theta;
    std::vector<Map::single_landmark_s> sensor_range_landmarks;
    double weight = 1.0;
    
    //std::cout << "particle number-----" <<num_particles << "----i------" << i<< std::endl;
    
    vector<int> associations;
    vector<double> sense_x;
    vector<double> sense_y;
    
    // List of Measurements in sensor range
    for (unsigned int j=0;j<map_landmarks.landmark_list.size();j++)
    {
      double xl = map_landmarks.landmark_list[j].x_f;
      double yl = map_landmarks.landmark_list[j].y_f;
      if (dist(xp, yp, xl,yl) < sensor_range)
      {
         sensor_range_landmarks.push_back(map_landmarks.landmark_list[j]);
      }        
    }
    
    
    for (unsigned int j=0;j<observations.size();j++)
    {
      
      double xc = observations[j].x;
      double yc = observations[j].y;
      double min_dist = std::numeric_limits<double>::max();
      
      int lm_id;
      double lm_x, lm_y;
      

      // transform car observation to map coordinate
      double xm = xp + xc * cos(thetap) - yc * sin(thetap);
      double ym = yp + xc * sin(thetap) + yc * cos(thetap);
      //std::cout << "observation number-----" <<observations.size() << "----j------" << j << std::endl;
      
      for (unsigned int k=0;k<sensor_range_landmarks.size();k++)
      {
         double xl = sensor_range_landmarks[k].x_f;
         double yl = sensor_range_landmarks[k].y_f;
        
         double distance = dist(xm, ym, xl, yl);
         
         if (min_dist > distance)
         {
           min_dist = distance; 
           lm_id = sensor_range_landmarks[k].id_i;
            lm_x = xl;
            lm_y = yl;
         }
         //std::cout << "sensor_range_landmarks-----" <<sensor_range_landmarks.size() << "---k------" << k<< std::endl;
      
      }
      associations.push_back(lm_id);
      sense_x.push_back(lm_x);
      sense_y.push_back(lm_y);
      
      weight *= multiv_prob(std_landmark[0], std_landmark[1], xm, ym, lm_x, lm_y);
     
  }
  // set the associations for particle "i"
  SetAssociations(particles[i], associations, sense_x, sense_y);
  // set weight for particle "i"
  particles[i].weight = weight;
  //std::cout << "particle filter updateWeights " << particles[i].weight << "......." << std::endl;
  // add particle weight to weight list
  weights.push_back(weight);

 }
}
                                     
                                     


void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   
  // discrete distribution based on particle weights
  std::discrete_distribution<> dis_dist(weights.begin(), weights.end());

  //random number engine
  std::default_random_engine gen;

  // resample vector of particles
  vector<Particle> sample_particles (particles.size());
  for(int i = 0; i < num_particles; i++){
     sample_particles[i] = particles[dis_dist(gen)];
     std::cout << "resample filter updateWeights " << sample_particles[i].weight << "......." << std::endl;
    //std::cout << "particle filter updateWeights " << particles[i].weight << "......." << std::endl;
  }
  // update particles
  particles = sample_particles;
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