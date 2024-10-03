function interference_init

% Distributed control system: interference node
%
% Generates disturbing network traffic.

% Initialize TrueTime kernel
ttInitKernel('prioFP'); % nbrOfInputs, nbrOfOutputs, fixed priority

% Create sender task
starttime = 0;
period = 0.001;
prio = 1;
ttCreatePeriodicTask('interf_task', starttime, period, 'interfcode');

% Initialize network
ttCreateHandler('nw_handler', prio, 'msgRcvInterf');

