function [exectime, data] = actcode(seg, data)

switch seg
 case 1
  data.u = ttGetMsg;
  exectime = 0.045*rand;
 case 2
  ttAnalogOut(1, data.u)
  ttAnalogOut(2, data.u)
  ttAnalogOut(3, data.u)%四旋翼多加的
  ttAnalogOut(4, data.u)%四旋翼多加的
  ttAnalogOut(5, data.u)%四旋翼多加的
  ttAnalogOut(6, data.u)%四旋翼多加的
  exectime = -1; % finished
end
