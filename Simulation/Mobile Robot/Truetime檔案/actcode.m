function [exectime, data] = actcode(seg, data)

switch seg
 case 1
  data.u = ttGetMsg;
  exectime = 0.045*rand;% �o��D�n�]�w�����P�ʾ��һݪ��ɶ�����
 case 2
  ttAnalogOut(1, data.u(1));
  ttAnalogOut(2, data.u(2));
  exectime = -1; % finished
end
