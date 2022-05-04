import math

def main():

    # angle needed in degrees
    target_angle = 127.091194
    negative_target_angle = target_angle - 360

    # how close we want the angle to be, in degrees
    accuracy = 1

    earth_file = open("earth_positions_12Y.txt")
    target_file = open("jupiter_position_30.txt")

    day_count = 0

    write_file = open("earth_jupiter_angles.txt", 'w')

    for earth_line in earth_file:
        target_line = target_file.readline()

        # each of these lists contain timestep, x, y and z position, in that order
        earth_pos = earth_line.split()
        target_pos = target_line.split()

        # theta = arctan(Y/X)
        theta_earth = math.atan2(float(earth_pos[2]), float(earth_pos[1]))
        theta_target = math.atan2(float(target_pos[2]), float(target_pos[1]))

        theta = math.degrees(theta_target) - math.degrees(theta_earth)
        if target_angle - accuracy <= theta <= target_angle + accuracy or \
                negative_target_angle - accuracy <= theta <= negative_target_angle + accuracy:
            day_count += 1
            print("Day {} is a candidate, angle: {}".format(int(earth_pos[0]), theta))

        write_file.write("{}\t\t{}\n".format(int(earth_pos[0]), theta))

    write_file.close()
    print("Done. {} possible launch days found.".format(day_count))

    return

if __name__ == '__main__':
    main()
