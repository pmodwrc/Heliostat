# IO Module

This class handles the whole communication with the Heliostat via a CAN Bus from Kvaser. The class initializes the motors and encoders and starts sending signals back and fourth to update the `GUI` class. The function that waits for incoming messages and proceeds them is run in a single Thread so that it doesn't affect the other ongoing functions. On this way it is possible to send and receive messages from the CAN Bus at the same time. 

The commands that are given to the CAN Bus must be in Byte format. They are given as hexadecimal numbers. To find more information about the specific commands to move the motor and get the replies read the [manual](https://www.analog.com/media/en/dsp-documentation/software-manuals/TMCM-1060_TMCL_firmware_manual.pdf?isDownload=true) from the stepper motor. To get the specific CAN Bus communication protocoll for the rotatory encoders read the [manual](https://www.baumer.com/ch/de/p/30274). 

Additionally this class handles the information that the `Q4_com_module` gets from the 4Q sensor. This is better described on this [site](fourq.md).

In the folling part the [class](https://github.com/snowstorm26/Heliostat/blob/main/IO_module.py) itself is described in detail.

## Class documentation
:::IO_module.IO
