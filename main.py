import math
import sys

import nbody, calc_rel_velocity, horizons_api, calc_earth_angle as np

class Parameters:
    def __init__(self):
        self.N         = 0              # Number of particles
        self.tEnd      = 0              # time at which simulation ends
        self.dt        = 0              # timestep
        self.recordValues = 1           # position, velocity, acceleration are recorded at every recordValues'th timestep
        self.sail_area = 0              # surface area of the sail
        self.sail_mass = 0              # mass of the sail
        self.sail_vel_x = 0             # initial velocity of sail, in x
        self.sail_vel_y = 0             # initial velocity of sail, in y
        self.theta = 0                  # tilt angle of sail(in radians)
        self.target_name = ""           # target name
        self.target_mass = 0            # target mass
        self.target_velocity_x = 0      # target velocity in x
        self.target_velocity_y = 0      # target velocity in y
        self.target_x = 0               # target position in x
        self.target_y = 0               # target position in y
        self.target_distance =  0       # distance of target from Sun
        self.laser_distance = 0         # distance to which the laser will be used, for Starshot only
        self.science_radius = 0         # proximity to target in which the probe can do science, for Starshot only
        self.power = 0                  # power of laser
        self.starshot = False           # whether we are doing a starshot simulation or not
        self.output = True              # whether to print/write to file or not

def read_config(filename):
    file = open(filename, 'r')
    params = Parameters()

    file.readline()
    params.N = int(file.readline())
    file.readline()
    params.tEnd     =  float(file.readline().strip())
    file.readline()
    params.dt        = int(file.readline().strip())
    file.readline()
    params.recordValues = int(file.readline().strip())
    file.readline()
    params.sail_area = int(file.readline().strip())
    file.readline()
    params.sail_mass = float(file.readline())
    file.readline()
    params.sail_vel_x = float(file.readline())
    file.readline()
    params.sail_vel_y = float(file.readline())
    file.readline()
    params.theta = int(file.readline())
    file.readline()
    params.target_name = file.readline().strip()
    file.readline()
    params.target_mass = float(file.readline())
    file.readline()
    params.target_velocity_x = float(file.readline())
    file.readline()
    params.target_velocity_y = float(file.readline())
    file.readline()
    params.target_x = float(file.readline())
    file.readline()
    params.target_y = float(file.readline())
    file.readline()
    params.target_distance = float(file.readline())
    file.readline()
    params.laser_distance = float(file.readline())
    file.readline()
    params.science_radius = float(file.readline())
    file.readline()
    params.power = float(file.readline())

    file.close()
    return params


def main():
    if len(sys.argv) != 2:
        print("ERROR: config file name required")
        exit(1)

    # get the name of the config file
    config = sys.argv[1]

    # first, read config files to see what we are running
    params = read_config(config)

    # initialize parameters for the simulation
    if params.target_name == "Proxima Centauri":
        params.starshot = True

    c = 3e8

    # run simulation
    # if starshot
    if params.starshot:
        params.sail_vel_y = calc_rel_velocity.calc_vel(params.sail_mass, params.power/c, params.laser_distance) * c
        # supress output
        params.output = False
        # run, get target coordinates
        target_pos = nbody.run(params)
        params.target_x = target_pos[0]
        params.target_y = target_pos[1]
        # allow output
        params.output = True
        # run normally
        nbody.run(params)
    else:
        # run normally
        nbody.run(params)

    return

if __name__ == '__main__':
    main()
