function [exectime, data] = msgRcvActuator(~, data)

ttCreateJob('act_task')
disp('ERROR: sensor received a message');
exectime = -1;
