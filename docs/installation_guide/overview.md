# Installation guide

To install the Heliostat on a new PC connect the CAN Bus from the Heliostat via the D sub cable and the 4Q sensor via the serial port with the PC. Afterwards the the CAN Bus and the 4Q sensor should be visible in the device manager. Then install the appropriate driver from [Kvaser](https://www.kvaser.com/download/). The driver can also be found on the server of PMOD/WRC. Afterwards you have to get the executable file and the belonging ui files. Make sure the ui files are in a folder called `ui_files` that is in the same directory as the executable file. 

In case you want to start the Heliostat over the shell open the anaconda prompt and start the belonging environement with the following command.

'''
conda activate helio
'''

Afterwards start the python script main.py

Good luck and have fun with the Heliostat. 