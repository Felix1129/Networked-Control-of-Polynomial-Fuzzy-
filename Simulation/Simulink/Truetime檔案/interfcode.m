function [exectime, data] = interfcode(seg, data)

BWshare = 0.0; % Fraction of the network bandwidth occupied by this node
if (rand(1) < BWshare) 
  ttSendMsg(1, 1, 80);  % send 80 bits to myself
end
exectime = -1;

