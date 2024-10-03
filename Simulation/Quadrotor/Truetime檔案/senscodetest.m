function [exectime, data] = senscodetest(seg, data)

switch seg
 case 1
  data.x1 = ttAnalogIn(1);
  data.x2 = ttAnalogIn(2);
  data.x3 = ttAnalogIn(3);
  data.x4 = ttAnalogIn(4);%四旋翼多加的
  data.x5 = ttAnalogIn(5);%四旋翼多加的
  data.x6 = ttAnalogIn(6);%四旋翼多加的
  data.x = [data.x1;data.x2;data.x3;data.x4;data.x5;data.x6];
  exectime = 0.000001; % finished
 case 2
  ttSendMsg(2, data.x , 80); % Send message (80 bits) to node 2 (controller)
  exectime = 0.045*rand;
 case 3
  exectime = -1; % finished
end

