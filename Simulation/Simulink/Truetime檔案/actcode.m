function [exectime, data] = actcode(seg, data)

switch seg
 case 1
  data.u = ttGetMsg;
  exectime = 0.045*rand;% 這邊主要設定控制器到致動器所需的時間延遲
 case 2
  ttAnalogOut(1, data.u(1));
  ttAnalogOut(2, data.u(2));
  exectime = -1; % finished
end
