import requests

def main():
    parameters = {
        "format": "text",
        "COMMAND": "\'599\'", # 599 is Jupiter
        "OBJ_DATA": "NO",
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "VECTORS",
        "CENTER": "\'500@0\'", # Solar System Barycenter
        "START_TIME": "\'2022-09-22\'",
        "STOP_TIME": "\'2034-09-22\'",
        "STEP_SIZE": "\'1d\'",
        "VEC_TABLE": "\'2\'",
        "CSV_FORMAT": "YES"
    }

    response = requests.get("https://ssd.jpl.nasa.gov/api/horizons.api", params=parameters)

    print(response.status_code)
    print(response.reason)

    KM_TO_METER = 1000
    content = response.text.split("\n")
    processing_data = False
    file = open("horizons_jupiter.txt", 'w')
    for line in content:
        if line.startswith("$$SOE") or line.startswith("$$EOE"):
            processing_data = not processing_data
            continue
        if not processing_data:
            continue
        line = line.split(",")

        new_list = [float(x) * KM_TO_METER for x in [line[2], line[3]]]
        file.write("{:.3e}\t\t{:.3e}\n".format(new_list[0], new_list[1]))

    file.close()

    return

if __name__ == '__main__':
    main()
