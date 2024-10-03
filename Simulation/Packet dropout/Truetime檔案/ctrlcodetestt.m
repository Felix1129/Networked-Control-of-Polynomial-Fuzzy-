 function [exectime, data] = ctrlcodetest(seg, data)
global flag F11 F21 F12 F22

switch seg
    case 1
        y = ttGetMsg;       % Obtain sensor value
        x11=y(1);%%ex
        x22=y(2);%%ey
        x33=y(3);%%eth
        x44=y(4);%%vr
        x55=y(5);%%wr


        % % 前鑑部
        vr=x44;%%指定參考機器人之線速度
        wr=x55;%%指定參考機器人之角速度
        VR_max = 1;%%參考機器人之最大線速度
        WR_max = 1;%%參考機器人之最大角速度
        
        z1 = vr/ VR_max; %%第一前鑑部
        z2 = wr/abs(WR_max);%%第二前鑑部
        
        FOU1=0.1;
        FOU2=0.2;
        m = 0.4;
        n = 0.5;
        
        
        
        % 第一前鑑部歸屬度計算
     if z1 <= 0
        z11_up = 0;
        z12_up = 1;
     elseif z1 >=0 && z1 <=1
        z11_up = z1;
        z12_up = 1-z11_up;
     else
        z11_up = 1;
        z12_up = 0;
     end
        
       % 第一前鑑部第二歸屬度計算
    if z1 <= 0
        z11_low = 0;
        z12_low = 1;
    elseif z1 <= (0+FOU1)   
        z11_low = 0;
        z12_low = (1/(FOU1-1))*z1+1;
    elseif z1 <= (0+FOU1) && z1 <= (1-FOU1)
        z11_low = (1/(1-FOU1))*z1+((FOU1)/(FOU1-1));
        z12_low = (1/(FOU1-1))*z1+1;
    elseif z1 >= (1-FOU1) && z1 <= 1
        z11_low = (1/(1-FOU1))*z1+((FOU1)/(FOU1-1));
        z12_low = 0;
    else
        z11_low = 1;
        z12_low = 0;
    end

        % % 第二前鑑部歸屬度計算
     if z2 <=-1
        z21_up = 0;
        z22_up = 1;
     elseif z2 >=-1 && z2 <=1
        z21_up = 0.5*(z2+1);
        z22_up = 1-z21_up;
     else
        z21_up = 1;
        z22_up = 0;
     end
        
       % 第二前鑑部第二歸屬度計算
    if z2 <= -1
        z21_low = 0;
        z22_low = 1;
    elseif z2 <= (-1+FOU2)
        z21_low = 0;
        z22_low = (1/(FOU2-2))*z2+1;
    elseif z2 <= (0+FOU2) && z2 <= (1-FOU2)
        z21_low = (1/(2-FOU2))*z2+((FOU2)/(FOU2-2));
        z22_low = (1/(FOU2-2))*z2+1;
    elseif z2 >= (1-FOU2) && z2 <= 1
        z21_low = (1/(2-FOU2))*z2+((FOU2)/(FOU2-2));
        z22_low = 0;
    else
        z21_low = 1;
        z22_low = 0;
    end
        
    h1_up = z11_up*z21_up;
    h2_up = z12_up*z22_up;
    h1_low = z11_low*z21_low;
    h2_low = z12_low*z22_low;
    
    
    F11 = [0.64289	 -0.39278	-2.2401 ;
        0.10471	 0.18745	1.1284];
    F21 = [0.00021559	-0.00049962	-0.0026653	;
        0.00017444	-0.00048302	-0.0021647];
    F12 = [0.5684	-0.22918	-1.512	;
        0.09302	0.12389 	0.75833];
    F22 = [0.0000344	0.0000585	-0.00023651	;
        0.0000363	0.0000834   -0.00036147];
   

            fh1= m *  h1_low * F11(1,1) + n * h1_up * F11(1,1) + m *  h2_low * F21(1,1) + n * h2_up * F21(1,1);
            fh2= m *  h1_low * F11(1,2) + n * h1_up * F11(1,2) + m *  h2_low * F21(1,2) + n * h2_up * F21(1,2);
            fh3= m *  h1_low * F11(1,3) + n * h1_up * F11(1,3) + m *  h2_low * F21(1,3) + n * h2_up * F21(1,3);
            fh4= m *  h1_low * F11(2,1) + n * h1_up * F11(2,1) + m *  h2_low * F21(2,1) + n * h2_up * F21(2,1);
            fh5= m *  h1_low * F11(2,2) + n * h1_up * F11(2,2) + m *  h2_low * F21(2,2) + n * h2_up * F21(2,2);
            fh6= m *  h1_low * F11(2,3) + n * h1_up * F11(2,3) + m *  h2_low * F21(2,3) + n * h2_up * F21(2,3);
            
            
            %%控制器輸出v,w
            data.u1 =  (fh1 * (x11) + fh2 * (x22) + fh3 * (x33));%%控制器輸出v
            data.u2 =  (fh4 * (x11) + fh5 * (x22) + fh6 * (x33));%%控制器輸出w

            
        data.u = [data.u1;data.u2];
        %%%%
        ttAnalogOut(1,data.u(1));
        ttAnalogOut(2,data.u(2));
        ttAnalogOut(3,y(1));
        ttAnalogOut(4,y(2));
        ttAnalogOut(5,y(3));
        exectime = 0.000001; %very important so remember add exectime
    case 2
        ttSendMsg(2,data.u,80);% Send 80 bits to node 2 (actuator)
        exectime = -1; % finished
end
