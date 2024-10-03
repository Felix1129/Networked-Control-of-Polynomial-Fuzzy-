function [exectime, data] = ctrlcodetest(seg, data)

switch seg
 case 1
  y=ttGetMsg;       % Obtain sensor value
  data.u=1;
  ttAnalogOut(1,y(1));
  ttAnalogOut(2,y(2));
  ttAnalogOut(3,y(3));%�|���l�h�[��
  ttAnalogOut(4,y(4));%�|���l�h�[��
  ttAnalogOut(5,y(5));%�|���l�h�[��
  ttAnalogOut(6,y(6));%�|���l�h�[��
  exectime = 0.000001; %very important so remember add exectime
 case 2
  ttSendMsg(1, data.u, 80);    % Send 80 bits to node 1 (actuator)
  exectime = -1; % finished
end
