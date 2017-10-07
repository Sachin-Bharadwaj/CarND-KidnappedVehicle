/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

// set the number of particles, NOTE: this is tunable parameter
num_particles = 100;

// resize vector of particles and vector of weights
particles.resize(num_particles);
weights.resize(num_particles);

// Set the random number generator related settings
default_random_engine gen;

// create normal distribution for initial x,y,theta (i.e noise around GPS seed)
normal_distribution<double> dist_x(x, std[0]);
normal_distribution<double> dist_y(y, std[1]);
normal_distribution<double> dist_theta(theta, std[2]);

// Initialize all the particles using the normal distribution seeded by GPS estimate
for(int i=0; i < num_particles; i++){
  particles[i].id = i;
  particles[i].x = dist_x(gen);
  particles[i].y = dist_y(gen);
  particles[i].theta = dist_theta(gen);
  particles[i].weight = 1.0;
}

// set the initialized flag to true
is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  
  // Engine for later generation of particles
  default_random_engine gen;
  
  // Make distributions for adding noise
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  
  // Different equations based on if yaw_rate is zero or not
  for (int i = 0; i < num_particles; ++i) {
    
    if (abs(yaw_rate) != 0) {
      // Add measurements to particles
      particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
      particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
      particles[i].theta += yaw_rate * delta_t;
      
    } else {
      // Add measurements to particles
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
      // Theta will stay the same due to no yaw_rate
      
    }
    

    // Add noise to predictions
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
    } // end of for  

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

// define some constants to be used later while calculating the weights of the particles


}// end of function

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
const double cons1 = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
const double cons2 = 2 * std_landmark[0] * std_landmark[0];
const double cons3 = 2 * std_landmark[1] * std_landmark[1];

// loop over each particle
for (int i=0; i < num_particles; i++){

  // for finding the prod of multi variate Gaussian of each observation
  double prod = 1.0;

  // loop over each observation
  for (int j=0; j< observations.size(); j++) {
  
    // Transform the observation from vehicle cordinates to map cordinates from particle perspective
    double t_obs_x, t_obs_y;
    t_obs_x = observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta) + particles[i].x;
    t_obs_y = observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta) + particles[i].y;

    // Find the nearest landmark
    vector<Map::single_landmark_s> landmarks = map_landmarks.landmark_list;
    vector<double> landmark_obs_dist (landmarks.size());

    for (int k=0; k < landmarks.size(); k++) {
      // minimize the number of landmarks to search for in association using the sensor range
      double landmark_dist_int = sqrt(pow(particles[i].x - landmarks[k].x_f, 2) + pow(particles[i].y - landmarks[k].y_f, 2));
      // populate as is if distance is within sensor range
      if(landmark_dist_int <= sensor_range) {
        landmark_obs_dist[k] = sqrt(pow(t_obs_x - landmarks[k].x_f, 2) + pow(t_obs_y - landmarks[k].y_f, 2));
      }
      else { // put some large number
        landmark_obs_dist[k] = 1000000;
      }
    } // end of for k
    
    // perform association of the observation with nearest landmark neighbour
    int pos_min = distance(landmark_obs_dist.begin(),min_element(landmark_obs_dist.begin(),landmark_obs_dist.end()));
    float mu_x = landmarks[pos_min].x_f;
    float mu_y = landmarks[pos_min].y_f;
    
    // calculate the prod of gaussian
    prod *= cons1 * exp(-pow(t_obs_x - mu_x, 2)/cons2 - pow(t_obs_y - mu_y, 2)/cons3);
    } // end of for j
  
  // update particle weights
  particles[i].weight = prod;
  weights[i] = particles[i].weight;
   
  } // end of for i


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
// Vector for new particles
  vector<Particle> new_particles (num_particles);
  
  // Use discrete distribution to return particles by weight
  default_random_engine gen;
  for (int i = 0; i < num_particles; ++i) {
    discrete_distribution<int> index(weights.begin(), weights.end());
    new_particles[i] = particles[index(gen)];
    
  }
  
  // Replace old particles with the resampled particles
  particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
