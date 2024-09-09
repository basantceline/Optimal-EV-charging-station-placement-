Created by JY Wong, 26 October 2019
MATLAB program and Simulink model for the IEEE 33 bus system 

Reference:
Network reconfiguration in distribution
systems for loss reduction and load balancing
by Baran and Wu, 1989


* NOTES *

This model has been designed to complement the data provided by MATPOWER
(Version 7.0) in the case study named 'case33bw'

Voltage measurements are phase to phase 
Current measurements are per phase


For simple simulations, please use the simple model.

For more comprehensive control, please use the comprehensive model, which is complemented with the MATLAB m file.


* HOW TO USE THE SIMULINK MODEL: IEEE33BusTestSystem *


* Method 1: Running the model through the MATLAB program (in the same directory)
-> Uncomment the "CB Inputs (MATLAB)" section and comment out the "CB Inputs (Manual)" section
-> Set the circuit breaker status inputs accordingly
-> Run the MATLAB program (which will subsequently run the Simulink model) and acquire the output variables
 -- simoutV => Voltage measurements at each branch
 -- simoutI => Current measurements at each branch

* Method 2: Manually set the circuit breakers
-> Uncomment the "CB Inputs (Manual)" section and comment out the "CB Inputs (MATLAB)" section
-> Manually setting the circuit breaker status in the respective blocks
-> Place your scopes or other measuring blocks at points of interest in the system
-> Simulate the Simulink model and observe the outputs 



* FEEDBACK *

Hope this program benefits you in your work! 
Please comment below if you have any suggestions for rooms of improvements

