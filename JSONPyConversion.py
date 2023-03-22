import json
import data_dict
# FILE NAMES WILL CHANGE AS WE PROCEED. JUST FOR TESTING PURPOSES ATM


# can be used to convert Unity inputs to Dict values
def json_to_py():
    # Converts JSON file to a Dict type and creates a Py file with those keys/values
    with open('jsonEx.json') as json_file:
        data = json.load(json_file)
        with open('JsonToPy.py', 'w') as f:
            for key, value in data.items():
                f.write('%s = %s\n' % (key, value))


# (NOT FINISHED: NEEDS EXPECTED OUTPUTS)
# can be used to convert Py outputs into Json values for Unity
def py_to_json():
    # dummy dictionary to export whole file as a JSON file 1:1
    main_dict = {
        'meas_cons': data_dict.meas_cons,
        'det_cons': data_dict.det_cons,
        'lis_cons': data_dict.lis_cons,
        'given_cons': data_dict.given_cons
    }
    # creates translated file from py dict values
    with open('PyToJson.json', 'w') as f:
        json.dump(main_dict, f, indent=4)


if __name__ == '__main__':
    print("Json to Py: Open JsonToPy.py")
    json_to_py()
    print("Py to Json: Open PyToJson.json")
    py_to_json()
