## Sensors
### Rotatory encoders
Both Axes are equipped with 18 bit single turn rotary encoders. This would allow passive tracking within an accuracy below 0.01 degrees. This accuracy will be high enough for most radiometry experiments. The rotary encoders are programmed to automatically send the actual position through the CAN network. The rotary encoders are powered independently
from the motors with 24V DC.


### Optical sensors
In addition to the rotary encoders an active guiding system with a sun sensor is installed. This system allows a higher tracking accuracy. The optical sensor is aligned with the experiment in the laboratory and ensures the tracking stability directly where it is needed. The active tracking system is also independent from small misalignments in the mechanical arrangement. To date the signals from the optical sensors are acquired trough an RS232 connection. 

### Limit switches
As a safety measure each motion axis is equipped with emergency end switches. These emergency switches cut down the power supply if a certain angle is reached, so that the motors stop immediately. If an emergency switch is triggered, it needs to be bypassed in order to move the motor again. Additionally there is also a limit switch in the heliostat primary axis in the forward direction. This limit switch is connected to the motor driver. If the limit switch is reached, the motor driver automatically stops the movement in the forward direction. It is however still possible to move the motor backwards without bypassing the limit switch. This limit
switch prevents the user from arduous recovery manoeuvres with bypassing emergency switches, if the heliostat is accidentally left running at the end of the operational time