function [exectime, data] = senscodetest(seg, data)

switch seg
 case 1
  data.x1 = ttAnalogIn(1);% 將編號1的類比值讀入轉成離散值
  data.x2 = ttAnalogIn(2);% 將編號2的類比值讀入轉成離散值
  data.x3 = ttAnalogIn(3);% 將編號3的類比值讀入轉成離散值
  data.x4 = ttAnalogIn(4);% 將編號3的類比值讀入轉成離散值
  data.x5 = ttAnalogIn(5);% 將編號3的類比值讀入轉成離散值
  data.x = [data.x1;data.x2;data.x3;data.x4;data.x5];% 轉變成狀態向量的形式
  exectime = 0.000001; % finished
  case 2
  ttSendMsg(3, data.x , 80); % Send message (80 bits) to node 3 (controller)
  exectime = 0.045*rand;%%感測器至控制器時間延遲
  
 case 3
  exectime = -1; % finished
end

