import sys


def main():
    file = open("proximity.txt", 'r')
    min_proximity = sys.maxsize
    step = 0

    for line in file:
        split_line = line.split()
        distance = float(split_line[1])
        if min_proximity > distance:
            min_proximity = distance
            step = split_line[0]

    print("Closest proximity: {:.3e}\nTimestep: {}".format(min_proximity, step))
    return

if __name__ == '__main__':
    main()
