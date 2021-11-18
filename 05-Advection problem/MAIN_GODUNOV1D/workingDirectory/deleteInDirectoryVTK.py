import os

dir_name = "./"
myDir = os.listdir(dir_name)

for item in myDir:
    if item.endswith(".vtk"):
        os.remove(os.path.join(dir_name, item))