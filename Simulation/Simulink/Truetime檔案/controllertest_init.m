function controllertest_init(arg)

% Distributed control system: controller node
%
% Receives messages from the sensor node, computes control signal
% and sends it to the actuator node. Also contains a high-priority
% disturbing task.

% Initialize TrueTime kernel
ttInitKernel('prioFP'); % nbrOfInputs, nbrOfOutputs, fixed priority

% Controller parameters
h = 0.050;

% Create task data (local memory)
 data.u = 0.0;
% data.k1 = 1;
% data.k2 = 1;
% data.k3 = 1;

% Create controller task
deadline = h;
prio = 2;
ttCreateTask('tsfuzzy_task', deadline, 'ctrlcodetest', data);
ttAttachNetworkHandler('tsfuzzy_task')

% Optional disturbance task

% Initialize network
ttCreateHandler('nw_handler', prio, 'msgRcvCtrltest');
