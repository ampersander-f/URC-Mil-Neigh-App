'''
This file is meant to be a sort-of "main" method, where the simulation parameters are outlined and run. Simulations and analysis only need this file to be modified, ideally this should streamline the simulations across many parameter sets and make parallel computing easier.
'''

#imports
from analyzer import Analyzer
from PressureSensorRL import PressureSensorRL
from PressureSensor import PressureSensor
from Grapher import Grapher

import math
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi

#DASH location on supercomputing cluster
sys.path = sys.path + ['/home/cmliepold/Dans_GPU/md_engine-master7/build/python/build/lib.linux-x86_64-2.7']


#simulation parameters
temp = .9 			#reduced temp units
density = .93				#reduced density units

sigma = 52.98 			#angstroms
epsilon = 275.7 		#kT
mass = 758798.8 

# -----> Defining the simulation parameters
#reduced units parameter values
# sigma, epsilon, alpha, r_a, a_g, sig_g, r_g
param_list = [(1, 1, .11268, 1.1051, .96, 1.00547, .05153)]	#parameters after being reduced, used for RL
param_list = (epsilon, sigma, mass)	
box_params = (110.26 * sigma, 173.27*sigma, 110.26*sigma)						#size of the box


#setting up the simulation
command_list = [(100, 4000), (0, 100, 531250)]
#how to run the simulations themselves

#graphing method
g = Grapher("convolution")

#run the simulation
sim = SimulatorRL(box_params, density, temperature, 1, Vector(-.0166,0,0), param_list[0], command_list)

fn = "RL_caurs" #output name
sim.update_fn(fn)

#let's go, babey !
sim.run()



#----------> Analysis

g = Grapher("CAURS") #choose graphing style

file = fn + ".xyz"  #grab the output file from the simulation

#add the pressure sensors
ps = PressureSensorRL((50.1182*sigma, 60.1418*sigma, 76.6113*sigma, 96.6586*sigma), 2, param_list, g)
pressure_sensor_list = [ps]

#generate the analysis object instance 
ana = Analyzer( param_list, pressure_sensor_list, fn )

#time to run some science
f = open(file, "r")
number_of_particles = f.readline().strip()

#rip processor 
ana.perform_analysis(f, number_of_particles)

#and scene
f.close()


