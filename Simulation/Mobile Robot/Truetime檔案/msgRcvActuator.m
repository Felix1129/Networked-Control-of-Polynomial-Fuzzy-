function [exectime, data] = msgRcvActuator(seg, data)

ttCreateJob('act_task')
disp('ERROR: sensor received a message');
exectime = -1;
