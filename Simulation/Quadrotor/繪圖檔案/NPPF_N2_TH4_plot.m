 clear;                      % 清除Workdpsce之變量，釋放系統內存
clear all;                  % 清除simulink資料，否則會發生同名檔案重複使用
clc;                        % 清空Command Window之顯示內容
tic                         % 啟用秒錶計時器，記錄程式執行執行時間
close all
global m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 
Ix = 8.1*10^(-3);           % 宣告x軸之轉動慣量
Iy = 8.1*10^(-3);           % 宣告y軸之轉動慣量
Iz = 1.42*10^(-2);          % 宣告z軸之轉動慣量
b = 5.42*10^(-5);           % 宣告升力係數
d = 1.1*10^(-6);            % 宣告阻力係數
%12~12.5
t_start = 20;               % 宣告外部干擾發生初始時間
t_end = 21;
disturbance = -0.35;            % 宣告外部干擾大小
l = 0.24;                   % 宣告螺旋槳到四旋翼中心之距離
m = 1;                      % 宣告重量
g = 9.81;                   % 宣告重力加速度
simulation_time = '50';     % 宣告模擬總時間
model = 'NPPF_N2_TH4';          % 宣告Simulink檔案名稱

%%

th=1;




disturbance = -0.3;            % 宣告外部干擾大小


AL1 = [16.4411    7.5434];     % 宣告高度控制增益
AL2 = [ 8.9386    4.132];

P1 = [    5.7067        0   9.1851         0    1.6983         0; % 宣告位置控制增益
         0   66.7822         0   37.4722         0   21.4015];
P2 = [    2.4671         0    4.2063         0    0.7521         0;
         0   23.9894         0   14.0483         0    5.5111];

AT1 = [   14.5556   13.4824   -0.3469    2.6257    0.5422   -0.2918;% 宣告姿態控制增益
   12.7239   15.4217    1.0159    0.6329    2.5554    0.9328;
    1.0385    1.7125    2.2981   -0.1595    0.5284    1.6825];
AT2 = [   17.4330   16.5480   -0.2481    2.9632    0.8636   -0.2488;
   14.1968   17.0713    1.0816    0.8017    2.7648    0.9947;
    1.3238    2.0823    2.5468   -0.1637    0.6064    1.8650];

%%



sim(model)                                          % 進行動態系統模擬
set_param(model, 'StopTime', simulation_time)       % 設定模擬參數
%紀錄紅色資料，紅色線疊在外面
m1=x_1N(:,2);
m2=y_1N(:,2);
m3=z_1N(:,2);
m4=t_1N;
m5=ISEatt_1N;
m6=ISEpos_1N;
m7=ePsi_1N;
m8=eTheta_1N;
m9=ePhi_1N;
m10=ez_1N;
m11=ey_1N;
m12=ex_1N;
%%


for gg=1
    if th==1
        figure(1)                                                           % 創建圖形視窗
        hold on                                                             % 保留當前繪圖
        plot3(x_1N(:,1),y_1N(:,1),z_1N(:,1),'k--',...                         % 繪製圖形
            x_1N(:,2),y_1N(:,2),z_1N(:,2),'r','LineWidth',1.4)
        legend({'reference','Novel SOS N=2 ',...         % 於坐標區添加圖例
            'Existing SOS N=2 '},'FontSize',12)
        xlabel('x(m)')                                                      % 添加x軸標籤
        ylabel('y(m)')                                                      % 添加y軸標籤
        zlabel('z(m)')                                                      % 添加z軸標籤
        axis([-1.5  1.5  -1.5 1.5 0 6])
        set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
        box on                                                              % 顯示坐標邊框
        grid on                                                             % 顯示圖形內網格
   %     view(-15,20);                                                      % 設定圖形視角
        view(15,20);                                                      % 設定圖形視角
        
        
        figure(2)                                                           % 創建圖形視窗
        hold on                                                             % 保留當前繪圖
        plot(x_1N(:,1),y_1N(:,1),'k--',...                                    % 繪製圖形
            x_1N(:,2),y_1N(:,2),'r','LineWidth',1.4);
        legend({'reference','Novel SOS N=2 ',...         % 於坐標區添加圖例
            'Existing SOS N=2 '},'FontSize',12)
        xlabel('x(m)')                                                      % 添加x軸標籤
        ylabel('y(m)')                                                      % 添加y軸標籤
        axis([-2  3 -2 3])
        set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
        box on                                                              % 顯示坐標邊框
        grid on                                                             % 顯示圖形內網格
        
         
        figure(3)                                                           % 創建圖形視窗
        hold on                                                             % 保留當前繪圖
        plot(t_1N,ISEatt_1N,'r','LineWidth',1.4);                           % 繪製圖形
        legend({'Novel SOS N=2 ',...         % 於坐標區添加圖例
            'Existing SOS N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
        xlabel('Time(s)')                                                   % 添加x軸標籤
        ylabel('ISE of e_{att} (radian)')                                   % 添加y軸標籤
        axis([0  50 0 20])
        set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
        box on                                                              % 顯示坐標邊框
             
        
        figure(4)                                                           % 創建圖形視窗
        hold on                                                             % 保留當前繪圖
        plot(t_1N,ISEpos_1N,'r','LineWidth',1.4);                           % 繪製圖形
        legend({'Novel SOS N=2 ',...         % 於坐標區添加圖例
            'Existing SOS N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
        xlabel('Time(s)')                                                   % 添加x軸標籤
        ylabel('ISE of e_{pos} (m)')                                        % 添加y軸標籤
        axis([0  50 0 2])
        set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
        box on                                                              % 顯示坐標邊框
        
        figure(11)                                                          % 創建圖形視窗
        hold on                                                             % 保留當前繪圖
        plot(t_1N,ePsi_1N,'r','LineWidth',1.4);                             % 繪製圖形
        legend({'Novel SOS N=2 ',...         % 於坐標區添加圖例
            'Existing SOS N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
        xlabel('Time(s)')                                                   % 添加x軸標籤
        ylabel('e_\psi (radian)')                                           % 添加y軸標籤
        axis([0  50 -6 4])
        set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
        box on                                                              % 顯示坐標邊框
        
        
        figure(12)                                                          % 創建圖形視窗
        hold on                                                             % 保留當前繪圖
        plot(t_1N,eTheta_1N,'r','LineWidth',1.4);                           % 繪製圖形
        legend({'Novel SOS N=2 ',...         % 於坐標區添加圖例
            'Existing SOS N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
        xlabel('Time(s)')                                                   % 添加x軸標籤
        ylabel('e_\theta (radian)')                                         % 添加y軸標籤
        axis([0  50 -2 2])
        set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
        box on                                                              % 顯示坐標邊框
        
        
        figure(13)                                                          % 創建圖形視窗
        hold on                                                             % 保留當前繪圖
        plot(t_1N,ePhi_1N,'r','LineWidth',1.4);                             % 繪製圖形
        legend({'Novel SOS N=2 ',...         % 於坐標區添加圖例
            'Existing SOS N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
        xlabel('Time(s)')                                                   % 添加x軸標籤
        ylabel('e_\phi (radian)')                                           % 添加y軸標籤
        set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
        box on                                                              % 顯示坐標邊框
        
        
        figure(14)                                                          % 創建圖形視窗
        hold on                                                             % 保留當前繪圖
        plot(t_1N,ez_1N,'r','LineWidth',1.4);                               % 繪製圖形
        legend({'Novel SOS N=2 ',...         % 於坐標區添加圖例
            'Existing SOS N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
        xlabel('Time(s)')                                                   % 添加x軸標籤
        ylabel('e_z (m)')                                                   % 添加y軸標籤
                        axis([0  50 -0.1 0.2])
        set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
        box on                                                              % 顯示坐標邊框
        
        
        figure(15)                                                          % 創建圖形視窗
        hold on                                                             % 保留當前繪圖
        plot(t_1N,ey_1N,'r','LineWidth',1.4);                               % 繪製圖形
        legend({'Novel SOS N=2 ',...         % 於坐標區添加圖例
            'Existing SOS N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
        xlabel('Time(s)')                                                   % 添加x軸標籤
        ylabel('e_y (m)')                                                   % 添加y軸標籤
        axis([0  50 -0.2 0.4])
        set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
        box on                                                              % 顯示坐標邊框
        
        
        figure(16)                                                          % 創建圖形視窗
        hold on                                                             % 保留當前繪圖
        plot(t_1N,ex_1N,'r','LineWidth',1.4);                               % 繪製圖形
        legend({'Novel SOS N=2 ',...         % 於坐標區添加圖例
            'Existing SOS N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
        xlabel('Time(s)')                                                   % 添加x軸標籤
        ylabel('e_x (m)')                                                   % 添加y軸標籤
        axis([0  50 -2 2])
        set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
        box on                                                              % 顯示坐標邊框
        disturbance = -0.4;
       

    end
end



% %  th=2;
% % 
% % AL1 = [15.8231    7.2599];     % 宣告高度控制增益
% % AL2 = [   7.8574    3.6322];
% % 
% % P1 = [    5.1007         0    8.2097         0    1.5179         0; % 宣告位置控制增益
% %          0   59.6903         0   33.4928         0   19.1288];
% % P2 = [    2.2717         0    3.8732         0    0.6926         0
% %          0   22.0893         0   12.9356         0    5.0745];
% % 
% % AT1 = [   14.5556   13.4824   -0.3469    2.9758    0.6145   -0.3584;     % 宣告姿態控制增益
% %    12.7239   15.4217    1.0159    0.7173    2.8961    1.1457;
% %     1.0385    1.7125    2.2981   -0.1808    0.5988    2.0665];
% % AT2 = [   17.4330   16.5480   -0.2481    3.3583    0.9788   -0.3056;
% %    14.1968   17.0713    1.0816    0.9086    3.1334    1.2218;
% %     1.3238    2.0823    2.5468   -0.1855    0.6872    2.2907];
% % sim(model)                                          % 進行動態系統模擬
% % set_param(model, 'StopTime', simulation_time)       % 設定模擬參數
% % 
% % for gg=1
% %     if th==2
% %         figure(1)                                                           % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot3( x_1N(:,2),y_1N(:,2),z_1N(:,2),'b','LineWidth',1.4)% 繪製圖形
% %          legend({'reference','Novel IT2 N=2 ',...         % 於坐標區添加圖例
% %             'Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)
% %         xlabel('x(m)')                                                      % 添加x軸標籤
% %         ylabel('y(m)')                                                      % 添加y軸標籤
% %         zlabel('z(m)')                                                      % 添加z軸標籤
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         grid on                                                             % 顯示圖形內網格
% %         view(-15,30);                                                      % 設定圖形視角
% %         
% %         
% %         figure(2)                                                           % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot( x_1N(:,2),y_1N(:,2),'b','LineWidth',1.4);% 繪製圖形
% %          legend({'reference','Novel IT2 N=2 ',...         % 於坐標區添加圖例
% %             'Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)
% %         xlabel('x(m)')                                                      % 添加x軸標籤
% %         ylabel('y(m)')                                                      % 添加y軸標籤
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         grid on                                                             % 顯示圖形內網格
% %         
% %         
% %         figure(3)                                                           % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ISEatt_1N,'b','LineWidth',1.4);                           % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('ISE of e_{att} (radian)')                                   % 添加y軸標籤
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(4)                                                           % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ISEpos_1N,'b','LineWidth',1.4);                           % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('ISE of e_{pos} (m)')                                        % 添加y軸標籤
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(11)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ePsi_1N,'b','LineWidth',1.4);                             % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_\psi (radian)')                                           % 添加y軸標籤
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(12)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,eTheta_1N,'b','LineWidth',1.4);                           % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_\theta (radian)')                                         % 添加y軸標籤
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(13)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ePhi_1N,'b','LineWidth',1.4);                             % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_\phi (radian)')                                           % 添加y軸標籤
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(14)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ez_1N,'b','LineWidth',1.4);                               % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_z (m)')                                                   % 添加y軸標籤
% %                         axis([0  50 -0.1 0.2])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(15)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ey_1N,'b','LineWidth',1.4);                               % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_y (m)')                                                   % 添加y軸標籤
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(16)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ex_1N,'b','LineWidth',1.4);                               % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_x (m)')                                                   % 添加y軸標籤
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% % 
% %         disturbance = -0.5;
% %         
% % 
% %     end
% % end
% % 
% % th=3;
% % 
% % AL1 = [   13.1035    6.0121];    % 宣告高度控制增益
% % AL2 = [  6.8482    3.1657];
% % 
% % P1 = [    4.3937         0    7.0717         0    1.3075         0;     % 宣告位置控制增益
% %          0   51.4164         0   28.8502         0   16.4773];
% % P2 = [    1.9053         0    3.2485         0    0.5809         0;
% %          0   18.5265         0   10.8492         0    4.2561];
% % 
% % AT1 = [   14.5556   13.4824   -0.3469    3.5260    0.7281   -0.4409;  % 宣告姿態控制增益
% %    12.7239   15.4217    1.0159    0.8499    3.4315    1.4093;
% %     1.0385    1.7125    2.2981   -0.2142    0.7095    2.5420];
% % AT2 = [   17.4330   16.5480   -0.2481    3.9792    1.1597   -0.3759
% %    14.1968   17.0713    1.0816    1.0765    3.7127    1.5029
% %     1.3238    2.0823    2.5468   -0.2198    0.8143    2.8178]
% % 
% % sim(model)                                          % 進行動態系統模擬
% % set_param(model, 'StopTime', simulation_time)       % 設定模擬參數
% % 
% % for gg=1
% %     if th==3
% %         figure(1)                                                           % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot3( x_1N(:,2),y_1N(:,2),z_1N(:,2),'g','LineWidth',1.4)% 繪製圖形
% %         legend({'reference','Novel IT2 N=2 ',...         % 於坐標區添加圖例
% %             'Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)
% %         xlabel('x(m)')                                                      % 添加x軸標籤
% %         ylabel('y(m)')                                                      % 添加y軸標籤
% %         zlabel('z(m)')                                                      % 添加z軸標籤
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         grid on                                                             % 顯示圖形內網格
% %         view(-15,30);                                                      % 設定圖形視角
% %         
% %         
% %         figure(2)                                                           % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot( x_1N(:,2),y_1N(:,2),'g','LineWidth',1.4);% 繪製圖形
% %          legend({'reference','Novel IT2 N=2 ',...         % 於坐標區添加圖例
% %             'Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)
% %         xlabel('x(m)')                                                      % 添加x軸標籤
% %         ylabel('y(m)')                                                      % 添加y軸標籤
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         grid on                                                             % 顯示圖形內網格
% %    
% %         
% %         figure(3)                                                           % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ISEatt_1N,'g','LineWidth',1.4);                           % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('ISE of e_{att} (radian)')                                   % 添加y軸標籤
% %         axis([0  50 0 20])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %              
% %         
% %         figure(4)                                                           % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ISEpos_1N,'g','LineWidth',1.4);                           % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('ISE of e_{pos} (m)')                                        % 添加y軸標籤
% %         axis([0  50 0 10])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(11)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ePsi_1N,'g','LineWidth',1.4);                             % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_\psi (radian)')                                           % 添加y軸標籤
% %         axis([0  50 -6 4])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(12)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,eTheta_1N,'g','LineWidth',1.4);                           % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_\theta (radian)')                                         % 添加y軸標籤
% %         axis([0  50 -2 2])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(13)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ePhi_1N,'g','LineWidth',1.4);                             % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_\phi (radian)')                                           % 添加y軸標籤
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(14)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ez_1N,'g','LineWidth',1.4);                               % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_z (m)')                                                   % 添加y軸標籤
% %                         axis([0  50 -0.1 0.2])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(15)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ey_1N,'g','LineWidth',1.4);                               % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_y (m)')                                                   % 添加y軸標籤
% %         axis([0  50 -0.2 0.4])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(16)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(t_1N,ex_1N,'g','LineWidth',1.4);                               % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_x (m)')                                                   % 添加y軸標籤
% %         axis([0  50 -2 2])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %        
% %         th=4;
% %     end
% % end
% % 
% % 
% % for gg=1
% %     if th==4
% %         figure(1)                                                           % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot3( m1,m2,m3,'r','LineWidth',1.4)% 繪製圖形
% %         legend({'reference','Novel IT2 N=2 ',...         % 於坐標區添加圖例
% %             'Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)
% %         xlabel('x(m)')                                                      % 添加x軸標籤
% %         ylabel('y(m)')                                                      % 添加y軸標籤
% %         zlabel('z(m)')                                                      % 添加z軸標籤
% %                 axis([-1.5  1.5  -1.5 1.5 0 6])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         grid on                                                             % 顯示圖形內網格
% %         view(15,20);                                                       % 設定圖形視角
% %         
% %         
% %         figure(2)                                                           % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot( m1,m2,'r','LineWidth',1.4);% 繪製圖形
% %          legend({'reference','Novel IT2 N=2 ',...         % 於坐標區添加圖例
% %             'Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)
% %         xlabel('x(m)')                                                      % 添加x軸標籤
% %         ylabel('y(m)')                                                      % 添加y軸標籤
% %                 axis([-1.5  1.5 -1.5 1.5])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         grid on                                                             % 顯示圖形內網格
% %    
% %         
% %         figure(3)                                                           % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(m4,m5,'r','LineWidth',1.4);                           % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('ISE of e_{att} (radian)')                                   % 添加y軸標籤
% %         axis([0  50 0 1])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %              
% %         
% %         figure(4)                                                           % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(m4,m6,'r','LineWidth',1.4);                           % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('ISE of e_{pos} (m)')                                        % 添加y軸標籤
% %         axis([0  50 0 0.5])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %     
% %         
% %         figure(11)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(m4,m7,'r','LineWidth',1.4);                             % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_\psi (radian)')                                           % 添加y軸標籤
% %         axis([0  50 -1.1 0.5])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(12)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(m4,m8,'r','LineWidth',1.4);                           % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_\theta (radian)')                                         % 添加y軸標籤
% %         axis([0  50 -0.05 0.15])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(13)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(m4,m9,'r','LineWidth',1.4);                             % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_\phi (radian)')                                           % 添加y軸標籤
% %              axis([0  50 -0.02 0.1])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(14)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(m4,m10,'r','LineWidth',1.4);                               % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_z (m)')                                                   % 添加y軸標籤
% %                         axis([0  50 -0.2 0.3])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(15)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(m4,m11,'r','LineWidth',1.4);                               % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_y (m)')                                                   % 添加y軸標籤
% %         axis([0  50 -0.2 0.4])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% %         
% %         
% %         figure(16)                                                          % 創建圖形視窗
% %         hold on                                                             % 保留當前繪圖
% %         plot(m4,m12,'r','LineWidth',1.4);                               % 繪製圖形
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % 於坐標區添加圖例
% %         xlabel('Time(s)')                                                   % 添加x軸標籤
% %         ylabel('e_x (m)')                                                   % 添加y軸標籤
% %         axis([0  50 -0.3 0.5])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % 設定圖形屬性參數
% %         box on                                                              % 顯示坐標邊框
% % 
% %     end
% % end

%}
toc                                                                 % 顯示程式總執行時間