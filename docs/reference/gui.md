# GUI Module

In this module includes all classes that handle and controll the different plots and the GUI. The most important class is the `HeliostatUI` class.
This is the main class of the controll software and connects the frontend to the backend. This class handels the main window from which all functions and plots are controlled. Beside the initialization of the buttons and the variables it also has the function update state which can be seen as the main function.

This function is called every second and first of all asks for the current position of the mirrors and then switches the motors on depending on the current mode that the operator wants to execute. It calls the corresponding functions which in the end all calculate the speed for the two motors. Afterwards this data is sent to the motor and written into a log file via the `Logger.log()` function from the `Logger` class. Furthermore it updates the various plots which will be explained later. To prevent the motor from running into a position that could damage the cables or mirrors a limit check is done in the end. In case the motor reaches one of these positions it is turned of immediately.

In the following parts the different functions from the Heliostat are  explained. To get further information visit the [code](https://github.com/snowstorm26/Heliostat/blob/main/GUI.py).

## Class documentation
::: GUI.HeliostatUI
    








