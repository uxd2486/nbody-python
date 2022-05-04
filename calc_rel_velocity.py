import math

def calc_vel(sail_mass, force, x_final):
    x = 0
    t = 0
    v = 0
    acc_0 = force / sail_mass
    dt = 0.1 # in seconds
    c = 3e8

    while x < x_final:
        # calculates velocity as a multiple of c
        v = ((acc_0 * t) / c) / (math.sqrt(1 + (((acc_0 * t) / c) ** 2)))
        x += (v*c*dt)
        t += dt

    return v

def main():
    file = open("velocities.txt", 'w')
    file.write("Distance\tRelativistic Velocity\n")
    # in order, Moon-Earth, laser, Mars-Earth, Earth-Sun, Jupiter-Sun, all in meters
    test_distances = [3.844e8, 2.186e9, 1.85e10, 5.46e10, 1.5e11, 7.786e11]
    sail_mass = 1e-3  # 1 gram
    force = 2e11 / 3e8   # power = 10^11 Watts, force = 2*power / c

    for distance in test_distances:
        print("Distance: {:.3e}".format(distance))
        relative_vel = calc_vel(sail_mass, force, distance)
        print("Velocity: {:.2}c\n".format(relative_vel))
        file.write("{:.3e}\t{:.2}c\n".format(distance, relative_vel))

    print("Done.")
    file.close()
    return

if __name__ == '__main__':
    main()
