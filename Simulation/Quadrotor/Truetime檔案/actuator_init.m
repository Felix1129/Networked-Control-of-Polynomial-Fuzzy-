function actuator_init

% Distributed control system: actuator node
%
% Receives messages from the controller and actuates 
% the plant.

% Initialize TrueTime kernel
ttInitKernel('prioFP'); % nbrOfInputs, nbrOfOutputs, fixed priority

% Create sensor task
data.x1 =ttAnalogIn(1);
data.x2 =ttAnalogIn(2);
data.x3 =ttAnalogIn(3);
data.x4 =ttAnalogIn(4);%四旋翼多加的
data.x5 =ttAnalogIn(5);%四旋翼多加的
data.x6 =ttAnalogIn(6);%四旋翼多加的
starttime = 0;
period = 0.050;
ttCreatePeriodicTask('sens_task', starttime, period, 'senscodetest', data);

% Create actuator task
deadline = 100;
prio = 1;
ttCreateTask('act_task', deadline, prio, 'actcode');
ttAttachNetworkHandler('act_task')

% Initialize network
ttCreateHandler('nw_handler', prio, 'msgRcvActuator');
