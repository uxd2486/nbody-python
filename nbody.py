import numpy as np
import math
import time

"""
Create Your Own N-body Simulation (With Python)
Philip Mocz (2020) Princeton Univeristy, @PMocz

Simulate orbits of stars interacting due to gravity
Code calculates pairwise forces according to Newton's Law of Gravity
"""

def getGravAcc( pos, mass, G, softening ):
	"""
    Calculate the acceleration on each particle due to Newton's Law 
	pos  is an N x 3 matrix of positions
	mass is an N x 1 vector of masses
	G is Newton's Gravitational constant
	softening is the softening length
	a is N x 3 matrix of accelerations
	"""
	# positions r = [x,y,z] for all particles
	x = pos[:,0:1]
	y = pos[:,1:2]
	z = pos[:,2:3]

	# matrix that stores all pairwise particle separations: r_j - r_i
	dx = x.T - x
	dy = y.T - y
	dz = z.T - z

	# matrix that stores 1/r^3 for all particle pairwise particle separations 
	inv_r3 = (dx**2 + dy**2 + dz**2 + softening**2)
	inv_r3[inv_r3>0] = inv_r3[inv_r3>0]**(-1.5)

	ax = G * (dx * inv_r3) @ mass
	ay = G * (dy * inv_r3) @ mass
	az = G * (dz * inv_r3) @ mass
	
	# pack together the acceleration components
	a = np.hstack((ax,ay,az))

	return a

def getAcc(pos, mass, G, softening, sail_area, theta):
	"""
    Calculate the total acceleration on each particle, including acceleration
    due to photon force which is only experienced by the solar sail
	pos  is an N x 3 matrix of positions
	mass is an N x 1 vector of masses
	G is Newton's Gravitational constant
	softening is the softening length
	sail_area is the total surface area of the solar sail
	acc is N x 3 matrix of accelerations
	"""
	# get acceleration due to gravity from provided code
	acc = getGravAcc(pos, mass, G, softening)

	# get the sail parameters
	sail_pos = pos[1]
	sail_mass = mass[1][0]
	sail_grav_acc = acc[1]

	# tilt angle
	phi = np.arctan2(sail_pos[1],sail_pos[0])
	totalAngle = theta + phi

	effectiveArea = sail_area * np.cos(theta) * np.cos(theta)

	# calculate (A*10^17) / m
	num = (2 * 1e17 * effectiveArea) / sail_mass

	# calculate 1/r^2
	inv_r2 = ( sail_pos[0]**2 + sail_pos[1]**2 + sail_pos[2]**2 + softening**2 ) ** (-1.0)

	# calculate acceleration due to photon force
	ax = num * inv_r2 * np.cos(totalAngle)
	ay = num * inv_r2 * np.sin(totalAngle)
	az = 0
	phot_acc = np.hstack((ax,ay,az))

	# add accelerations together
	sail_acc = sail_grav_acc + phot_acc

	# set the correct acceleration for the sail
	acc[1] = sail_acc

	return acc

	
def getEnergy( pos, vel, mass, G ):
	"""
	Get kinetic energy (KE) and potential energy (PE) of simulation
	pos is N x 3 matrix of positions
	vel is N x 3 matrix of velocities
	mass is an N x 1 vector of masses
	G is Newton's Gravitational constant
	KE is the kinetic energy of the system
	PE is the potential energy of the system
	"""
	# Kinetic Energy:
	KE = 0.5 * np.sum(np.sum( mass * vel**2 ))

	# Potential Energy:

	# positions r = [x,y,z] for all particles
	x = pos[:,0:1]
	y = pos[:,1:2]
	z = pos[:,2:3]

	# matrix that stores all pairwise particle separations: r_j - r_i
	dx = x.T - x
	dy = y.T - y
	dz = z.T - z

	# matrix that stores 1/r for all particle pairwise particle separations 
	inv_r = np.sqrt(dx**2 + dy**2 + dz**2)
	inv_r[inv_r>0] = 1.0/inv_r[inv_r>0]

	# sum over upper triangle, to count each interaction only once
	PE = G * np.sum(np.sum(np.triu(-(mass*mass.T)*inv_r,1)))
	
	return KE, PE

def run(params):
	""" N-body simulation """

	# Simulation parameters
	t         = 0              # current time of the simulation
	G         = 6.67e-11       # Newton's Gravitational Constant
	event_horizon = 4*6.597e8  # radius from the sun beyond which our simulation goes haywire
	softening = 0.1            # softening length

	# to check if we are within the science radius of our target
	doing_science = False

	# to store start and end time of when we are close enough to do science
	science_start = 0
	science_end = 0

	# store the arrival timestep and position
	arrival_step = 0
	arrival_pos = 0

	# check if we have arrived yet
	arrived = False

	# Generate Initial Conditions
	mass = np.ones((params.N,1))    # initialize masses as 1
	pos  = np.zeros((params.N,3))   # initialize position vectors as 0
	vel  = np.zeros((params.N,3))   # initialize velocity vectors as 0

	# Set initial positions and velocities for all objects
	# For Sun
	mass[0][0] = 1.9891e30

	# For solar sail, initial position from Horizons at 24 Sep 2022
	mass[1][0] = params.sail_mass
	pos[1][0] = 1.487339502126168e11
	pos[1][1] = 1.789363230737144e9
	vel[1][0] = params.sail_vel_x
	vel[1][1] = params.sail_vel_y

	# target
	mass[2][0] = params.target_mass
	pos[2][0] = params.target_x
	pos[2][1] = params.target_y

	# Convert to Center-of-Mass frame
	vel -= np.mean(mass * vel,0) / np.mean(mass)

	# calculate initial gravitational accelerations
	acc = getAcc(pos, mass, G, softening, params.sail_area, params.theta)

	# calculate initial energy of system
	KE, PE  = getEnergy( pos, vel, mass, G )

	# number of timesteps
	Nt = int(np.ceil(params.tEnd/params.dt))

	# save energies, particle orbits for plotting trails
	pos_save = np.zeros((params.N,3,Nt+1))
	pos_save[:,:,0] = pos
	KE_save = np.zeros(Nt+1)
	KE_save[0] = KE
	PE_save = np.zeros(Nt+1)
	PE_save[0] = PE

	# store position/acceleration of sail for every iteration
	pos_file = open("sail_positions.txt".format(math.degrees(params.theta)), 'w')

	# write the timestep and initial total energy to the file
	if params.output:
		print("Beginning simulation.")
	initial_total_energy = KE+PE

	# record starting time
	start_time = time.time()

	# Simulation Main Loop
	for i in range(Nt):
		if params.output:
			print("Loop ",i)

		if i % params.recordValues == 0 and params.output:
			pos_file.write("{}\t\t{:.3e}\t\t{:.3e}\t\t{:.3e}\n".format(i,pos[1][0], pos[1][1], pos[1][2]))

		# (1/2) kick
		vel += acc * params.dt/2.0

		# drift
		pos += vel * params.dt

		# update accelerations
		acc = getAcc(pos, mass, G, softening, params.sail_area, params.theta)

		# (1/2) kick
		vel += acc * params.dt/2.0

		# update time
		t += params.dt

		# if we're too close to the Sun, stop the simulation
		distance_to_sun = math.sqrt((pos[1][0]**2) + (pos[1][1]**2))
		if distance_to_sun < event_horizon:
			break

		# check if the sail is at target orbit yet
		if not arrived:
			if distance_to_sun > params.target_distance:
				arrival_pos = [pos[1][0], pos[1][1]]
				arrival_step = i
				arrived = True

		# see how close to the target we are
		distance = math.sqrt((pos[2][0] - pos[1][0]) ** 2 + (pos[2][1] - pos[1][1]) ** 2)

		science_proximity =  distance < params.science_radius

		if not doing_science:
			# check if we can start doing science
			if science_proximity:
				doing_science = True
				science_start = i
		else:
			# check if we have to stop doing science
			if not science_proximity:
				doing_science = False
				science_end = i

		# get energy of system
		KE, PE  = getEnergy( pos, vel, mass, G )

		# save energies, positions for plotting trail
		pos_save[:,:,i+1] = pos
		KE_save[i+1] = KE
		PE_save[i+1] = PE

	# record ending time
	end_time = time.time()

	# get the final total energy
	final_total_energy = KE+PE

	# simulation is finished
	pos_file.close()
	if params.output:
		print("Simulation finished!")
		print("{} is the timestep at which sail reaches the target".format(arrival_step))
		print("Target position: {:.5e}\t{:.5e}".format(arrival_pos[0], arrival_pos[1]))
		print("Start: {}\nEnd: {}\nTime: {}".format(science_start, science_end, science_end-science_start))

	if params.output:
	# simulation.txt stores timestep, initial and final total energy and time taken
		energy_file = open("simulation_NFP.txt".format(params.dt), 'a')
		energy_file.write("Theta: {}\N{DEGREE SIGN}\n".format(math.degrees(params.theta)))
		energy_file.write("Timestep: {}\nInitial Total Energy: {:.5e} J\n".format(params.dt, initial_total_energy))
		energy_file.write("Final Total Energy: {:.5e} J\n".format(final_total_energy))
		energy_file.write("Change in Energy: {:.7e} J\n".format(final_total_energy-initial_total_energy))
		energy_file.write("Time taken: {:.5} seconds\n\n".format(end_time-start_time))
		energy_file.close()

	if params.starshot and not params.output:
		return arrival_pos

	return 0
