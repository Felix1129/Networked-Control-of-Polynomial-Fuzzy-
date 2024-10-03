 function [exectime, data] = ctrlcodetest(seg, data)
global flag

switch seg
    case 1
        y = ttGetMsg;       % Obtain sensor value
        x11=y(1);%%ex
        x22=y(2);%%ey
        x33=y(3);%%eth
        x44=y(4);%%vr
        x55=y(5);%%wr
        % %%SOS
        %{
        vr=x44;%%指定參考機器人之線速度
        wr=x55;%%指定參考機器人之角速度
        
        VR_max = 1;%%參考機器人之最大線速度
        WR_max = 1;%%參考機器人之最大角速度
        
        z1 = vr/ VR_max; %%第一前鑑部
        z2 = wr/abs(WR_max);%%第二前鑑部
        
        % %   第一前鑑部歸屬度計算
        if z1<=0
            vr_1=1;
            vr_2=0;
        elseif z1>=0 && z1<=1
            vr_1=1-z1;
            vr_2=1-vr_1;
        else
            vr_1=0;
            vr_2=1;
        end
        % %   第二前鑑部歸屬度計算
        if z2<=-1
            wr_1=1;
            wr_2=0;
            wr_3=0;
        elseif z2<=0 && z2>-1
            wr_2=z2+1;
            wr_1=1-wr_2;
            wr_3=0;
        else
            wr_2=1-z2;
            wr_3=1-wr_2;
            wr_1=0;
        end
        
        h1=vr_1*wr_1;
        h2=vr_1*wr_2;
        h3=vr_1*wr_3;
        h4=vr_2*wr_1;
        h5=vr_2*wr_2;
        h6=vr_2*wr_3;
        
        
        
        %{
B=gain();
for i=1:3
    for j=1:2
        F1(j,i)=B(j,i);
        F2(j,i)=B(j+2,i);
        F3(j,i)=B(j+4,i);
        F4(j,i)=B(j+6,i);
        F5(j,i)=B(j+8,i);
        F6(j,i)=B(j+10,i);
    end
end
        %}
        
        
        % % %D=1
        
        F1 =[-0.1621   -0.1350   -0.0606;
            0   -0.0473   -0.0647];
        F2 =[-1.2539   -0.3771   -0.3205;
            0   -0.8996   -1.3470];
        F3 =[ -1.1516    0.1265    0.7584;
            0.0001   -0.9359   -0.6818];
        F4 =[ -0.9248   -2.0755   -2.1529;
            -0.0086   -2.7000   -2.0093];
        F5 =[ -0.3240    1.6831   -0.1350;
            -0.0005   -1.0209   -1.0867];
        F6 =[-1.1217    1.3977    1.3595;
            0.0075   -1.2494   -2.2590];
        
        
        fh1=h1*F1(1,1)+h2*F2(1,1)+h3*F3(1,1)+h4*F4(1,1)+h5*F5(1,1)+h6*F6(1,1);
        fh2=h1*F1(1,2)+h2*F2(1,2)+h3*F3(1,2)+h4*F4(1,2)+h5*F5(1,2)+h6*F6(1,2);
        fh3=h1*F1(1,3)+h2*F2(1,3)+h3*F3(1,3)+h4*F4(1,3)+h5*F5(1,3)+h6*F6(1,3);
        fh4=h1*F1(2,1)+h2*F2(2,1)+h3*F3(2,1)+h4*F4(2,1)+h5*F5(2,1)+h6*F6(2,1);
        fh5=h1*F1(2,2)+h2*F2(2,2)+h3*F3(2,2)+h4*F4(2,2)+h5*F5(2,2)+h6*F6(2,2);
        fh6=h1*F1(2,3)+h2*F2(2,3)+h3*F3(2,3)+h4*F4(2,3)+h5*F5(2,3)+h6*F6(2,3);
        % % %控制器輸出v,w
        data.u1=-(fh1*x11+fh2*x22+fh3*x33);%%控制器輸出v
        data.u2=-(fh4*x11+fh5*x22+fh6*x33);%%控制器輸出w
        %}
        % % % % % % %%新二型多項式
        vr=x44;%%指定參考機器人之線速度
        wr=x55;%%指定參考機器人之角速度
        VR_max = 1;%%參考機器人之最大線速度
        WR_max = 1;%%參考機器人之最大角速度
        %%%%前鑑部
        z1 = vr/ VR_max; %%第一前鑑部
        z2 = wr/abs(WR_max);%%第二前鑑部
        
        v_up = z1;
        v_low = 1-z1;
        
        w_up = (z2+1)/2;
        w_low = (-z2+1)/2;
        % 第一前鑑部歸屬度計算
        if v_up >= 1
            v_up = 1;
        end
        if v_up <= 0
            v_up = 0;
        end
        if  v_low >= 1
            v_low = 0;
        end
        if  v_low <= 0
            v_low = 1;
        end
        %      % 第二前鑑部歸屬度計算
        if  w_up >= 1
            w_up = 1;
        end
        if  w_up <=-1
            w_up = 0;
        end
        if  w_low >= 1
            w_low = 0;
        end
        if  w_low <= -1
            w_low = 1;
        end
        h1=(w_up*2+v_up)/3;
        h2=(w_low*2+v_low)/3;
        
        if flag == 1
            F11 =[-2.3400   -2.6   -1.6404;
                -1.3   -9.1000   -5.000];
            F21 =[-2.3400    0.3120    0.1690;
                0.1950   -8.7750   -7.200];
        end
        if flag == 2
            F11 =[-2.3400   -2.6   -1.6404;
                -1.3   -9.1000   -5.500];
            F21 =[-2.3400    0.3120    0.1690;
                0.1950   -8.7750   -7.920];
        end
        if flag == 3
            F11 =[-2.3400   -2.6   -1.6404;
                -1.3   -9.1000   -6.125];
            F21 =[-2.3400    0.3120    0.1690;
                0.1950   -8.7750   -8.820];
        end
%             
%         G11 =[-1.8000   -2   -1.6404;
%             -1   -7.0000   -2.5000];
%         G21 =[-1.8000    0.2400    0.1300;
%             0.1500   -6.7500   -3.6000];
%         
%         

        
        
        % %D=1
        fh1=h1*F11(1,1)+h2*F21(1,1);
        fh2=h1*F11(1,2)+h2*F21(1,2);
        fh3=h1*F11(1,3)+h2*F21(1,3);
        fh4=h1*F11(2,1)+h2*F21(2,1);
        fh5=h1*F11(2,2)+h2*F21(2,2);
        fh6=h1*F11(2,3)+h2*F21(2,3);
        % %%控制器輸出v,w
        data.u1=-(fh1*x11+fh2*x22+fh3*x33);%%控制器輸出v
        data.u2=-(fh4*x11+fh5*x22+fh6*x33);%%控制器輸出w
        % % %%%傳輸控制訊號
        
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
