clc;
clear all;

name = 'wpt';
new_system(name);
open_system(name);

% solver definitions: default
add_block('powerlib/powergui',[name, '/powergui']);
set_param([name, '/powergui'],'SimulationMode','Discrete');
set_param([name, '/powergui'],'SampleTime','1e-6');
set_param('wpt', 'MaxStep', '1');

Transmitter(name, Hierarchy([0,0,200,200]), 0);
Receiver(name, Hierarchy([200,0,1000,200]), 0);