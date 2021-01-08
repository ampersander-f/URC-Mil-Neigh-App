'''
This is a defined area inside the simulation that mimics the pressure sensors ("Wilhelmy Plates") in the experiments, so as to replicate the experimental conditions as closely as possible. it also considerably reduces the number of nearest neighbors i have to calculate :3
'''

import math
from Grapher import Grapher
import numpy as np

class PressureSensor:

	'''Creates the sensor itself, only defining its location and bounds.
	Parameters: 
		bounds - tuple of (x1, x2, y1, y2)
		frame - float of how large around the frame you can calculate potentials
		unreduced_params - (epsilon, sigma, mass) before they are converted to be one, used to calculate real values for the isotherms generated
		g - Grapher object
	Returns:
		null
	'''
	def __init__(self, bounds, frame, unreduced_params, g):
		self.bounds = bounds
		self.frame_bounds = tuple([i + frame for i in bounds])
		self.ur_eps = unreduced_params[0]
		self.ur_sig = unreduced_params[1]
		self.ur_mass = unreduced_params[2]

		self.plate_area = self.calc_plate_area(bounds)

		#holder lists
		self.areas = []				#store the areas of the simulation at a timestep
		self.pressures = []			#store the pressures of the sensor at a timestep

		self.grapher = g
		self.z_bins = []


	'''Calculates the area of the plate, in non-reduced units
	Parameters:
		bounds - boundary of the plate in (x1, x2, y1, y2) tuple form
	Returns:
		the area, in real units
	'''
	def calc_plate_area(self, bounds):
		return ( abs(bounds[0] - bounds[1]) * abs(bounds[2] - bounds[3]) ) #* self.ur_sig * self.ur_sig * 1e-9


	'''Used to find all the particles in a given region.
	Parameters: 
		area_to_search - list of parsed particle positions from which the points are pulled from
		bounds - what area is searched, given as a tuple of (x1, x2, y1, y2)
	Returns: 
		list with all of the particles within the area, each as a tuple of (x, y, z)
	'''
	def find_particles_in(self, area_to_search, bounds):
		particles_in_region = []				#where to store the particles to return

		x1 = bounds[0]
		x2 = bounds[1]
		y1 = bounds[2]
		y2 = bounds[3]
		for particle in area_to_search:
			x = particle[0]
			y = particle[1]

			#check if in bounds
			if (x >= x1 and x <= x2 and y >= y1 and y <= y2):
				particles_in_region.append(particle)


		return particles_in_region


	'''Calculates the force, determined by the derivative of the LJ potential function and the parameters passed into this pressure sensor instance. Must be modified in all subclasses that use a different inter-particle interaction potential to correctly calculate the force.
	Parameters: 
		r: distance (float)
	Returns: 
		the force (float)
	'''
	def force_calc(self, r):
		sigma = 1

		r_m = 61.0			#bottom of well, in Angstroms

		#force_const = 24 * self.ur_eps * sigma_6

		if(r != float(0) ):
			#power = r ** -6
			#force = ( (2 * sigma_6 * r ** -1 * power * power) - (r ** -1 * power) ) * force_const
			return 12.0 * self.ur_eps / r * ( (r_m/r)**12 - (r_m/r)**6 )

		return 0.0


	'''Calculates the pressure felt on the sensor, and returns the value divided by the area of the sensor itself.
	Parameters:
		nn_dict - dictionary of particles and nearest (and next-nearest) neighbors that gets iterated over
		area_of_simulation - the area of the simulation at the time step, as a float in non-reduced units
	Returns:
		Null (updates object variables)
	'''
	def calc_pressure(self, nn_dict, area_of_simulation):
		l = [i for i in nn_dict.iterkeys()]
		plate_particles = self.find_particles_in(l, self.bounds)

		pressure_sum = [0,0]

		for particle in plate_particles:
			neighbors = nn_dict[particle]
			px = particle[0]
			py = particle[1]
			pz = particle[2]

			for n in neighbors:
				r = math.sqrt( ((n[0] - px) * (n[0] - px)) + ((n[1] - py) * (n[1] - py))  + ((n[2] - pz) * (n[2] - pz)) )
				rx = (n[0] - px)
				ry = (n[1] - py)
				force = self.force_calc(r)

				d = np.sqrt(rx *rx + ry * ry)
				forcex = force * (rx / d)
				forcey = force * (ry / d)

				pressure_sum[0] += forcex * rx
				pressure_sum[1] += forcey * ry

		pressure_sum = [i / (2. * self.plate_area) for i in pressure_sum]
		self.pressures.append(pressure_sum)
		self.areas.append(area_of_simulation)

	'''
	Calculates the potential function based on the pressures and areas from above (see paper for formulae); uses same inputs as above
	'''
	def calc_pot(self, nn_dict, area_of_simulation):
		l = [i for i in nn_dict.iterkeys()]
		plate_particles = self.find_particles_in(l, self.bounds)

		pot_sum = [0,0]
		p_sum = [0]

		for particle in plate_particles:
			neighbors = nn_dict[particle]
			px = particle[0]
			py = particle[1]
			pz = particle[2]

			for n in neighbors:
				rx = (n[0] - px)
				ry = (n[1] - py)
				potx = self.potential_calc(rx)			#this is only in the subclass?
				poty = self.potential_calc(ry)
				pot_sum[0] += potx
				pot_sum[1] += poty

		self.pressures.append(pot_sum)
		self.areas.append(area_of_simulation)
		



	'''Calls the grapher object.
	Parameters;
		fn - the file name of the graph for a single simulation
	Returns:
		null
	'''
	def generate_graph(self, fn):	

		self.areas = self.areas[10:]
		self.pressures = self.pressures[10:]

		self.grapher.update_area(self.areas)
		self.grapher.update_pressure(self.pressures)
		self.grapher.graph("Potential1", fn)






