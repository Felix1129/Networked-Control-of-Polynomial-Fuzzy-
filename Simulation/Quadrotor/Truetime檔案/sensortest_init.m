function sensortest_init

% Distributed control system: sensor node
%
% Samples the plant periodically and sends the samples to the 
% controller node.

% Initialize TrueTime kernel
ttInitKernel(3, 0, 'prioFP'); % nbrOfInputs, nbrOfOutputs, fixed priority

% Create sensor task
data.x1 =ttAnalogIn(1);
data.x2 =ttAnalogIn(2);
data.x3 =ttAnalogIn(3);
data.x4 =ttAnalogIn(4);%四旋翼多加的
data.x5 =ttAnalogIn(5);%四旋翼多加的
data.x6 =ttAnalogIn(6);%四旋翼多加的

offset = 0;
period = 0.050;
prio = 1;
ttCreatePeriodicTask('sens_task', offset, period, prio, 'senscodetest', data);

% Initialize network
ttCreateInterruptHandler('nw_handler', prio, 'msgRcvSensor');
ttInitNetwork(4, 'nw_handler'); % node #4 in the network
