function [exectime, data] = senscodetest(seg, data)

switch seg
 case 1
  data.x1 = ttAnalogIn(1);% �N�s��1�������Ū�J�ন������
  data.x2 = ttAnalogIn(2);% �N�s��2�������Ū�J�ন������
  data.x3 = ttAnalogIn(3);% �N�s��3�������Ū�J�ন������
  data.x4 = ttAnalogIn(4);% �N�s��3�������Ū�J�ন������
  data.x5 = ttAnalogIn(5);% �N�s��3�������Ū�J�ন������
  data.x = [data.x1;data.x2;data.x3;data.x4;data.x5];% ���ܦ����A�V�q���Φ�
  exectime = 0.000001; % finished
  case 2
  ttSendMsg(3, data.x , 80); % Send message (80 bits) to node 3 (controller)
  exectime = 0.045*rand;%%�P�����ܱ���ɶ�����
  
 case 3
  exectime = -1; % finished
end

