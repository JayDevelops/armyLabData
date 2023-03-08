import os
import docker
from pyodm import Node, types
import time
import zipfile

# TODO:
#  - fix pathing to work with relevant folder(s) (no full path)
#  - add cases if more status errors occur


def main():
    path = r'C:..\ODMScripts\ODMImages'
    dir_list = os.listdir(path)
    path_list = []
    last_percentage = -1

    # each image path in ODMImages is added to a list
    for x in dir_list:
        with open(os.path.join(path, x), 'rb') as fh:
            path_list.append(fh.name)

    # docker starts and ports to the ODM server
    client = docker.from_env()
    container = client.containers.run("opendronemap/nodeodm", ports={"3000/tcp": 3000}, detach=True)

    # wait for container to fully connect
    time.sleep(4)

    # send all the images in ODMImages to the ODM server for 3D model creation
    n = Node('localhost', 3000)
    task = n.create_task(path_list)

    # display current 3D model progress
    while task.info().status == types.TaskStatus.RUNNING:
        # if the % is not the same as previously displayed, display it
        if task.info().progress != last_percentage:
            print(str(task.info().progress) + "% Complete")
            last_percentage = task.info().progress

    # download results from the server into a zipped 'results' folder and extract it
    zip_path = task.download_zip("results")
    with zipfile.ZipFile(zip_path, "r") as zip_h:
        zip_h.extractall("results")

    # delete the zipped file and stop the container
    os.remove(zip_path)
    container.stop()


if __name__ == '__main__':
    main()
