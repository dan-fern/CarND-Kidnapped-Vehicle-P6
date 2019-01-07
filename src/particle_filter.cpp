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


// random engine across multiple and various method calls
static default_random_engine gen;


void ParticleFilter::init(
			double x,
			double y,
			double theta,
			double std[ ] )
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// TODO: Set the number of particles
	num_particles = 101;

	// define normal sensor noise distributions
	normal_distribution<double> noise_x( 0, std[ 0 ] );
	normal_distribution<double> noise_y( 0, std[ 1 ] );
	normal_distribution<double> noise_theta( 0, std[ 2 ] );

	// initialize particles and add noise
	for( int i = 0; i < num_particles; i++ )
	{
		Particle p;
		p.id = i;
		p.x = x + noise_x( gen );
		p.y = y + noise_y( gen );
		p.theta = theta + noise_theta( gen );
		p.weight = 1.0;

		particles.push_back(p);
	}

	is_initialized = true;
}


void ParticleFilter::prediction(
			double delta_t,
			double std_pos[ ],
			double velocity,
			double yaw_rate )
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// define normal sensor noise distributions
    normal_distribution<double> noise_x( 0, std_pos[ 0 ] );
    normal_distribution<double> noise_y( 0, std_pos[ 1 ] );
    normal_distribution<double> noise_theta( 0, std_pos[ 2 ] );

    for( int i = 0; i < num_particles; i++ )
	{

	// calculate new state
	if( fabs( yaw_rate ) < 0.000005 )
	{
		particles[ i ].x += velocity * delta_t * cos( particles[ i ].theta );
		particles[ i ].y += velocity * delta_t * sin( particles[ i ].theta );
	}
	else
	{
		particles[ i ].x += velocity / yaw_rate *
				( sin( particles[ i ].theta + yaw_rate * delta_t ) -
				sin( particles[ i ].theta ) );
		particles[ i ].y += velocity / yaw_rate *
				( cos( particles[ i ].theta ) -
				cos( particles[ i ].theta + yaw_rate * delta_t ) );
		particles[ i ].theta += yaw_rate * delta_t;
	}

	// add noise
	particles[ i ].x += noise_x( gen );
	particles[ i ].y += noise_y( gen );
	particles[ i ].theta += noise_theta( gen );

	}

}


void ParticleFilter::dataAssociation(
			std::vector<LandmarkObs> predicted,
			std::vector<LandmarkObs>& observations )
{
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	for( unsigned int i = 0; i < observations.size( ); i++ )
	{

		// get current observation
		LandmarkObs obs = observations[ i ];

		// initialize minimum distance to maximum possible
		double minimum_distance = numeric_limits<double>::max( );

		// initialize associated landmark id from map placeholder
		int map_id = -1;

		for( unsigned int j = 0; j < predicted.size( ); j++ )
		{
			// get current prediction
			LandmarkObs pred = predicted[ j ];

			// get distance between observed and predicted landmarks
			double current_distance = dist( obs.x, obs.y, pred.x, pred.y );

			// find predicted landmark nearest the current observed landmark
			if ( current_distance < minimum_distance )
			{
				minimum_distance = current_distance;
				map_id = pred.id;
			}
		}

		// set observation id to the nearest predicted landmark id
		observations[ i ].id = map_id;
	}
}


void ParticleFilter::updateWeights(
			double sensor_range,
			double std_landmark[ ],
			const std::vector<LandmarkObs> &observations,
			const Map &map_landmarks )
{
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

	for( int i = 0; i < num_particles; i++ )
	{

		// get the particle x, y coordinates and theta
		double p_x = particles[ i ].x;
		double p_y = particles[ i ].y;
		double p_theta = particles[ i ].theta;

		// vector to hold landmark locations predicted to be within sensor range
		vector<LandmarkObs> predictions;

		// for each landmark...
		for( unsigned int j = 0; j < map_landmarks.landmark_list.size( ); j++ )
		{

			// get landmark id and x, y coordinates
			float l_x = map_landmarks.landmark_list[ j ].x_f;
			float l_y = map_landmarks.landmark_list[ j ].y_f;
			int l_id = map_landmarks.landmark_list[ j ].id_i;

			// only consider landmarks within sensor range of the particle
			if( fabs( l_x - p_x ) <= sensor_range &&
					fabs( l_y - p_y ) <= sensor_range )
			{

			  	// add prediction to vector
			  	predictions.push_back( LandmarkObs{ l_id, l_x, l_y } );
			}
		}

		// create copy of the list of observations
		vector<LandmarkObs> transformed_observations;

		// transform from vehicle coordinates to map coordinates
		for( unsigned int j = 0; j < observations.size( ); j++ )
		{

			double t_x = cos( p_theta ) * observations[ j ].x -
					sin( p_theta ) * observations[ j ].y + p_x;

			double t_y = sin( p_theta ) * observations[ j ].x +
					cos( p_theta ) * observations[ j ].y + p_y;

			transformed_observations.push_back(
					LandmarkObs{ observations[ j ].id, t_x, t_y } );
		}

		// call dataAssociation for the predictions and transformed observations
		dataAssociation( predictions, transformed_observations );

		// reinitilialize weight
		particles[ i ].weight = 1.0;

		for( unsigned int j = 0; j < transformed_observations.size( ); j++)
		{

			// initialize observation and associated prediction coordinates
			double obs_x, obs_y, pred_x, pred_y;
			obs_x = transformed_observations[ j ].x;
			obs_y = transformed_observations[ j ].y;

			int associated_prediction = transformed_observations[ j ].id;

			// get x, y coordinates of the associated prediction
			for( unsigned int k = 0; k < predictions.size( ); k++ )
			{
			  	if( predictions[ k ].id == associated_prediction )
				{
			    	pred_x = predictions[ k ].x;
			    	pred_y = predictions[ k ].y;
			  	}
			}

			// calculate observation weight with multivariate Gaussian
			double s_x = std_landmark[ 0 ];
			double s_y = std_landmark[ 1 ];

			double obs_w = ( 1 / ( 2 * M_PI * s_x * s_y ) ) *
					exp( -( pow( pred_x - obs_x, 2 ) / ( 2 * pow( s_x, 2 ) ) +
					( pow( pred_y - obs_y, 2 ) / ( 2 * pow( s_y, 2 ) ) ) ) );

			// product of this observation weight with total observations weight
			particles[ i ].weight *= obs_w;
		}
    }
}


void ParticleFilter::resample( )
{
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    vector<Particle> new_particles;

    // get all of the current weights
    vector<double> weights;

	for( int i = 0; i < num_particles; i++ )
	{
      	weights.push_back( particles[ i ].weight );
	}

	// generate random starting index for resampling wheel
	uniform_int_distribution<int> intDist( 0, num_particles - 1 );
	auto index = intDist( gen );

	// get max weight
	double max_weight = *max_element( weights.begin( ), weights.end( ) );

	// uniform random distribution [0.0, max_weight)
	uniform_real_distribution<double> realDist( 0.0, max_weight );

	double beta = 0.0;

	// resampling loop
	for( int i = 0; i < num_particles; i++ )
	{
	  	beta += realDist( gen ) * 2.0;
	  	while( beta > weights[ index ] )
		{
			beta -= weights[ index ];
			index = ( index + 1 ) % num_particles;
	  	}

	  	new_particles.push_back( particles[ index ] );
	}

	particles = new_particles;
}


void ParticleFilter::SetAssociations(
			Particle& particle,
			const std::vector<int>& associations,
            const std::vector<double>& sense_x,
			const std::vector<double>& sense_y )
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}


string ParticleFilter::getAssociations( Particle best )
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}


string ParticleFilter::getSenseX( Particle best )
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}


string ParticleFilter::getSenseY( Particle best )
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
