function [exectime, data] = ctrlcodetest(seg, data)

switch seg
 case 1
  y=ttGetMsg;       % Obtain sensor value
  data.u=1;
  ttAnalogOut(1,y(1));
  ttAnalogOut(2,y(2));
  ttAnalogOut(3,y(3));%四旋翼多加的
  ttAnalogOut(4,y(4));%四旋翼多加的
  ttAnalogOut(5,y(5));%四旋翼多加的
  ttAnalogOut(6,y(6));%四旋翼多加的
  exectime = 0.000001; %very important so remember add exectime
 case 2
  ttSendMsg(1, data.u, 80);    % Send 80 bits to node 1 (actuator)
  exectime = -1; % finished
end
