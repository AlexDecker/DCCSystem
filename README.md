# DCCSystem
DCC - Do not Carry Cables system. It's a distributed algorithm for controlling a wireless power transfet network. It's implemented using sdWPTSimulator (https://github.com/AlexDecker/sdWPT-Simulator)

## Usage
Open the MATLAB environment and type 
```Matlab
start
```
Example of generating new data from config1.m, dummie\_tx.m, the empty RX application, with parameters defined by params1.m, for all mobData from mobData1 and considering 10 repetitions
```Matlab
generateData("dummie", "params1", "configs1", [1,2,3,4,5,6,7,8,9,10], 10)
```
Exemple of plotting comparative charts torward the files bla.mat, ble.mat and bli.mat (using cell array)
```Matlab
plotComparativeData({'bla.mat','ble.mat','bli.mat'})
```

## Implemented systems/protocols/applications
 * dummie: uses random voltage vectors

## About the directory structure and files

 * "Configs" foder contains functions that defines some configurations regarding the environment and devices

 * The "main" folder has two files. "generateData" defines a function responsible for getting the configurations from the specified files and for running the simulation using the specified mobility data and number of repetitions. "simulate", in turn, is a general simulation script (see sdWPT-Simulator) that runs a single simulation using the specified applications, application params, hardware/environment configuration and mobility data.

 * The "mobility" folder contains the files with mobility data, grouped by the parameters that generated them probabilisticaly (using trackGenerator). The parameters are described on desc.txt (at each folder).

 * "optimizers" folder is dedicaded to the mathematical programming algorithms developed for solvig problems related to the operation of the recharging system.

 * "out" folder is dedicated to the output files of "generateData" function

 * "params" folder is dedicated to the functions that inform a set of parameters for the execution of the embedded applications.

 * "rxApplications"/"txApplications": see powerTXApplications and powerRXApplications at sdWPT-Simulator repository. They define the behavior of the devices and represent the embedded applications. The choice of which application to use on each device occurs in the call of the "simulate" (or "generateData") function. For the definition of new applications: 
 	* Copy/paste one of the already defined applications
	* Modify what you want
	* Modify the "appBuilder" function
	* Create a "param" file at "params/"

 * "utils":
 	* adjustAndSum: used by "generateData" for collecting some statistics about the data logged by the simulations
	* appBuilder: creates instances of embedded applications according the informed parameters and app name
	* evalSolution: evaluates quickly an informed solution in order to estimate the finishing time of the recharging process
	* plotComparativeData: compares the output data sets from "out"
