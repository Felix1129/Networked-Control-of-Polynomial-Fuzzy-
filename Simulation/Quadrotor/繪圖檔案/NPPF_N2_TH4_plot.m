 clear;                      % �M��Workdpsce���ܶq�A����t�Τ��s
clear all;                  % �M��simulink��ơA�_�h�|�o�ͦP�W�ɮ׭��ƨϥ�
clc;                        % �M��Command Window����ܤ��e
tic                         % �ҥά���p�ɾ��A�O���{���������ɶ�
close all
global m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 
Ix = 8.1*10^(-3);           % �ŧix�b����ʺD�q
Iy = 8.1*10^(-3);           % �ŧiy�b����ʺD�q
Iz = 1.42*10^(-2);          % �ŧiz�b����ʺD�q
b = 5.42*10^(-5);           % �ŧi�ɤO�Y��
d = 1.1*10^(-6);            % �ŧi���O�Y��
%12~12.5
t_start = 20;               % �ŧi�~���z�Z�o�ͪ�l�ɶ�
t_end = 21;
disturbance = -0.35;            % �ŧi�~���z�Z�j�p
l = 0.24;                   % �ŧi���ۼը�|���l���ߤ��Z��
m = 1;                      % �ŧi���q
g = 9.81;                   % �ŧi���O�[�t��
simulation_time = '50';     % �ŧi�����`�ɶ�
model = 'NPPF_N2_TH4';          % �ŧiSimulink�ɮצW��

%%

th=1;




disturbance = -0.3;            % �ŧi�~���z�Z�j�p


AL1 = [16.4411    7.5434];     % �ŧi���ױ���W�q
AL2 = [ 8.9386    4.132];

P1 = [    5.7067        0   9.1851         0    1.6983         0; % �ŧi��m����W�q
         0   66.7822         0   37.4722         0   21.4015];
P2 = [    2.4671         0    4.2063         0    0.7521         0;
         0   23.9894         0   14.0483         0    5.5111];

AT1 = [   14.5556   13.4824   -0.3469    2.6257    0.5422   -0.2918;% �ŧi���A����W�q
   12.7239   15.4217    1.0159    0.6329    2.5554    0.9328;
    1.0385    1.7125    2.2981   -0.1595    0.5284    1.6825];
AT2 = [   17.4330   16.5480   -0.2481    2.9632    0.8636   -0.2488;
   14.1968   17.0713    1.0816    0.8017    2.7648    0.9947;
    1.3238    2.0823    2.5468   -0.1637    0.6064    1.8650];

%%



sim(model)                                          % �i��ʺA�t�μ���
set_param(model, 'StopTime', simulation_time)       % �]�w�����Ѽ�
%���������ơA����u�|�b�~��
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
        figure(1)                                                           % �Ыعϧε���
        hold on                                                             % �O�d��eø��
        plot3(x_1N(:,1),y_1N(:,1),z_1N(:,1),'k--',...                         % ø�s�ϧ�
            x_1N(:,2),y_1N(:,2),z_1N(:,2),'r','LineWidth',1.4)
        legend({'reference','Novel SOS N=2 ',...         % �󧤼аϲK�[�Ϩ�
            'Existing SOS N=2 '},'FontSize',12)
        xlabel('x(m)')                                                      % �K�[x�b����
        ylabel('y(m)')                                                      % �K�[y�b����
        zlabel('z(m)')                                                      % �K�[z�b����
        axis([-1.5  1.5  -1.5 1.5 0 6])
        set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
        box on                                                              % ��ܧ������
        grid on                                                             % ��ܹϧΤ�����
   %     view(-15,20);                                                      % �]�w�ϧε���
        view(15,20);                                                      % �]�w�ϧε���
        
        
        figure(2)                                                           % �Ыعϧε���
        hold on                                                             % �O�d��eø��
        plot(x_1N(:,1),y_1N(:,1),'k--',...                                    % ø�s�ϧ�
            x_1N(:,2),y_1N(:,2),'r','LineWidth',1.4);
        legend({'reference','Novel SOS N=2 ',...         % �󧤼аϲK�[�Ϩ�
            'Existing SOS N=2 '},'FontSize',12)
        xlabel('x(m)')                                                      % �K�[x�b����
        ylabel('y(m)')                                                      % �K�[y�b����
        axis([-2  3 -2 3])
        set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
        box on                                                              % ��ܧ������
        grid on                                                             % ��ܹϧΤ�����
        
         
        figure(3)                                                           % �Ыعϧε���
        hold on                                                             % �O�d��eø��
        plot(t_1N,ISEatt_1N,'r','LineWidth',1.4);                           % ø�s�ϧ�
        legend({'Novel SOS N=2 ',...         % �󧤼аϲK�[�Ϩ�
            'Existing SOS N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
        xlabel('Time(s)')                                                   % �K�[x�b����
        ylabel('ISE of e_{att} (radian)')                                   % �K�[y�b����
        axis([0  50 0 20])
        set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
        box on                                                              % ��ܧ������
             
        
        figure(4)                                                           % �Ыعϧε���
        hold on                                                             % �O�d��eø��
        plot(t_1N,ISEpos_1N,'r','LineWidth',1.4);                           % ø�s�ϧ�
        legend({'Novel SOS N=2 ',...         % �󧤼аϲK�[�Ϩ�
            'Existing SOS N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
        xlabel('Time(s)')                                                   % �K�[x�b����
        ylabel('ISE of e_{pos} (m)')                                        % �K�[y�b����
        axis([0  50 0 2])
        set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
        box on                                                              % ��ܧ������
        
        figure(11)                                                          % �Ыعϧε���
        hold on                                                             % �O�d��eø��
        plot(t_1N,ePsi_1N,'r','LineWidth',1.4);                             % ø�s�ϧ�
        legend({'Novel SOS N=2 ',...         % �󧤼аϲK�[�Ϩ�
            'Existing SOS N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
        xlabel('Time(s)')                                                   % �K�[x�b����
        ylabel('e_\psi (radian)')                                           % �K�[y�b����
        axis([0  50 -6 4])
        set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
        box on                                                              % ��ܧ������
        
        
        figure(12)                                                          % �Ыعϧε���
        hold on                                                             % �O�d��eø��
        plot(t_1N,eTheta_1N,'r','LineWidth',1.4);                           % ø�s�ϧ�
        legend({'Novel SOS N=2 ',...         % �󧤼аϲK�[�Ϩ�
            'Existing SOS N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
        xlabel('Time(s)')                                                   % �K�[x�b����
        ylabel('e_\theta (radian)')                                         % �K�[y�b����
        axis([0  50 -2 2])
        set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
        box on                                                              % ��ܧ������
        
        
        figure(13)                                                          % �Ыعϧε���
        hold on                                                             % �O�d��eø��
        plot(t_1N,ePhi_1N,'r','LineWidth',1.4);                             % ø�s�ϧ�
        legend({'Novel SOS N=2 ',...         % �󧤼аϲK�[�Ϩ�
            'Existing SOS N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
        xlabel('Time(s)')                                                   % �K�[x�b����
        ylabel('e_\phi (radian)')                                           % �K�[y�b����
        set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
        box on                                                              % ��ܧ������
        
        
        figure(14)                                                          % �Ыعϧε���
        hold on                                                             % �O�d��eø��
        plot(t_1N,ez_1N,'r','LineWidth',1.4);                               % ø�s�ϧ�
        legend({'Novel SOS N=2 ',...         % �󧤼аϲK�[�Ϩ�
            'Existing SOS N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
        xlabel('Time(s)')                                                   % �K�[x�b����
        ylabel('e_z (m)')                                                   % �K�[y�b����
                        axis([0  50 -0.1 0.2])
        set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
        box on                                                              % ��ܧ������
        
        
        figure(15)                                                          % �Ыعϧε���
        hold on                                                             % �O�d��eø��
        plot(t_1N,ey_1N,'r','LineWidth',1.4);                               % ø�s�ϧ�
        legend({'Novel SOS N=2 ',...         % �󧤼аϲK�[�Ϩ�
            'Existing SOS N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
        xlabel('Time(s)')                                                   % �K�[x�b����
        ylabel('e_y (m)')                                                   % �K�[y�b����
        axis([0  50 -0.2 0.4])
        set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
        box on                                                              % ��ܧ������
        
        
        figure(16)                                                          % �Ыعϧε���
        hold on                                                             % �O�d��eø��
        plot(t_1N,ex_1N,'r','LineWidth',1.4);                               % ø�s�ϧ�
        legend({'Novel SOS N=2 ',...         % �󧤼аϲK�[�Ϩ�
            'Existing SOS N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
        xlabel('Time(s)')                                                   % �K�[x�b����
        ylabel('e_x (m)')                                                   % �K�[y�b����
        axis([0  50 -2 2])
        set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
        box on                                                              % ��ܧ������
        disturbance = -0.4;
       

    end
end



% %  th=2;
% % 
% % AL1 = [15.8231    7.2599];     % �ŧi���ױ���W�q
% % AL2 = [   7.8574    3.6322];
% % 
% % P1 = [    5.1007         0    8.2097         0    1.5179         0; % �ŧi��m����W�q
% %          0   59.6903         0   33.4928         0   19.1288];
% % P2 = [    2.2717         0    3.8732         0    0.6926         0
% %          0   22.0893         0   12.9356         0    5.0745];
% % 
% % AT1 = [   14.5556   13.4824   -0.3469    2.9758    0.6145   -0.3584;     % �ŧi���A����W�q
% %    12.7239   15.4217    1.0159    0.7173    2.8961    1.1457;
% %     1.0385    1.7125    2.2981   -0.1808    0.5988    2.0665];
% % AT2 = [   17.4330   16.5480   -0.2481    3.3583    0.9788   -0.3056;
% %    14.1968   17.0713    1.0816    0.9086    3.1334    1.2218;
% %     1.3238    2.0823    2.5468   -0.1855    0.6872    2.2907];
% % sim(model)                                          % �i��ʺA�t�μ���
% % set_param(model, 'StopTime', simulation_time)       % �]�w�����Ѽ�
% % 
% % for gg=1
% %     if th==2
% %         figure(1)                                                           % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot3( x_1N(:,2),y_1N(:,2),z_1N(:,2),'b','LineWidth',1.4)% ø�s�ϧ�
% %          legend({'reference','Novel IT2 N=2 ',...         % �󧤼аϲK�[�Ϩ�
% %             'Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)
% %         xlabel('x(m)')                                                      % �K�[x�b����
% %         ylabel('y(m)')                                                      % �K�[y�b����
% %         zlabel('z(m)')                                                      % �K�[z�b����
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         grid on                                                             % ��ܹϧΤ�����
% %         view(-15,30);                                                      % �]�w�ϧε���
% %         
% %         
% %         figure(2)                                                           % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot( x_1N(:,2),y_1N(:,2),'b','LineWidth',1.4);% ø�s�ϧ�
% %          legend({'reference','Novel IT2 N=2 ',...         % �󧤼аϲK�[�Ϩ�
% %             'Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)
% %         xlabel('x(m)')                                                      % �K�[x�b����
% %         ylabel('y(m)')                                                      % �K�[y�b����
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         grid on                                                             % ��ܹϧΤ�����
% %         
% %         
% %         figure(3)                                                           % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ISEatt_1N,'b','LineWidth',1.4);                           % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('ISE of e_{att} (radian)')                                   % �K�[y�b����
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(4)                                                           % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ISEpos_1N,'b','LineWidth',1.4);                           % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('ISE of e_{pos} (m)')                                        % �K�[y�b����
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(11)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ePsi_1N,'b','LineWidth',1.4);                             % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_\psi (radian)')                                           % �K�[y�b����
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(12)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,eTheta_1N,'b','LineWidth',1.4);                           % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_\theta (radian)')                                         % �K�[y�b����
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(13)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ePhi_1N,'b','LineWidth',1.4);                             % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_\phi (radian)')                                           % �K�[y�b����
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(14)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ez_1N,'b','LineWidth',1.4);                               % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_z (m)')                                                   % �K�[y�b����
% %                         axis([0  50 -0.1 0.2])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(15)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ey_1N,'b','LineWidth',1.4);                               % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_y (m)')                                                   % �K�[y�b����
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(16)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ex_1N,'b','LineWidth',1.4);                               % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_x (m)')                                                   % �K�[y�b����
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% % 
% %         disturbance = -0.5;
% %         
% % 
% %     end
% % end
% % 
% % th=3;
% % 
% % AL1 = [   13.1035    6.0121];    % �ŧi���ױ���W�q
% % AL2 = [  6.8482    3.1657];
% % 
% % P1 = [    4.3937         0    7.0717         0    1.3075         0;     % �ŧi��m����W�q
% %          0   51.4164         0   28.8502         0   16.4773];
% % P2 = [    1.9053         0    3.2485         0    0.5809         0;
% %          0   18.5265         0   10.8492         0    4.2561];
% % 
% % AT1 = [   14.5556   13.4824   -0.3469    3.5260    0.7281   -0.4409;  % �ŧi���A����W�q
% %    12.7239   15.4217    1.0159    0.8499    3.4315    1.4093;
% %     1.0385    1.7125    2.2981   -0.2142    0.7095    2.5420];
% % AT2 = [   17.4330   16.5480   -0.2481    3.9792    1.1597   -0.3759
% %    14.1968   17.0713    1.0816    1.0765    3.7127    1.5029
% %     1.3238    2.0823    2.5468   -0.2198    0.8143    2.8178]
% % 
% % sim(model)                                          % �i��ʺA�t�μ���
% % set_param(model, 'StopTime', simulation_time)       % �]�w�����Ѽ�
% % 
% % for gg=1
% %     if th==3
% %         figure(1)                                                           % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot3( x_1N(:,2),y_1N(:,2),z_1N(:,2),'g','LineWidth',1.4)% ø�s�ϧ�
% %         legend({'reference','Novel IT2 N=2 ',...         % �󧤼аϲK�[�Ϩ�
% %             'Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)
% %         xlabel('x(m)')                                                      % �K�[x�b����
% %         ylabel('y(m)')                                                      % �K�[y�b����
% %         zlabel('z(m)')                                                      % �K�[z�b����
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         grid on                                                             % ��ܹϧΤ�����
% %         view(-15,30);                                                      % �]�w�ϧε���
% %         
% %         
% %         figure(2)                                                           % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot( x_1N(:,2),y_1N(:,2),'g','LineWidth',1.4);% ø�s�ϧ�
% %          legend({'reference','Novel IT2 N=2 ',...         % �󧤼аϲK�[�Ϩ�
% %             'Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)
% %         xlabel('x(m)')                                                      % �K�[x�b����
% %         ylabel('y(m)')                                                      % �K�[y�b����
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         grid on                                                             % ��ܹϧΤ�����
% %    
% %         
% %         figure(3)                                                           % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ISEatt_1N,'g','LineWidth',1.4);                           % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('ISE of e_{att} (radian)')                                   % �K�[y�b����
% %         axis([0  50 0 20])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %              
% %         
% %         figure(4)                                                           % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ISEpos_1N,'g','LineWidth',1.4);                           % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('ISE of e_{pos} (m)')                                        % �K�[y�b����
% %         axis([0  50 0 10])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(11)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ePsi_1N,'g','LineWidth',1.4);                             % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_\psi (radian)')                                           % �K�[y�b����
% %         axis([0  50 -6 4])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(12)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,eTheta_1N,'g','LineWidth',1.4);                           % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_\theta (radian)')                                         % �K�[y�b����
% %         axis([0  50 -2 2])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(13)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ePhi_1N,'g','LineWidth',1.4);                             % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_\phi (radian)')                                           % �K�[y�b����
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(14)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ez_1N,'g','LineWidth',1.4);                               % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_z (m)')                                                   % �K�[y�b����
% %                         axis([0  50 -0.1 0.2])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(15)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ey_1N,'g','LineWidth',1.4);                               % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_y (m)')                                                   % �K�[y�b����
% %         axis([0  50 -0.2 0.4])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(16)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(t_1N,ex_1N,'g','LineWidth',1.4);                               % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_x (m)')                                                   % �K�[y�b����
% %         axis([0  50 -2 2])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %        
% %         th=4;
% %     end
% % end
% % 
% % 
% % for gg=1
% %     if th==4
% %         figure(1)                                                           % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot3( m1,m2,m3,'r','LineWidth',1.4)% ø�s�ϧ�
% %         legend({'reference','Novel IT2 N=2 ',...         % �󧤼аϲK�[�Ϩ�
% %             'Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)
% %         xlabel('x(m)')                                                      % �K�[x�b����
% %         ylabel('y(m)')                                                      % �K�[y�b����
% %         zlabel('z(m)')                                                      % �K�[z�b����
% %                 axis([-1.5  1.5  -1.5 1.5 0 6])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         grid on                                                             % ��ܹϧΤ�����
% %         view(15,20);                                                       % �]�w�ϧε���
% %         
% %         
% %         figure(2)                                                           % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot( m1,m2,'r','LineWidth',1.4);% ø�s�ϧ�
% %          legend({'reference','Novel IT2 N=2 ',...         % �󧤼аϲK�[�Ϩ�
% %             'Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)
% %         xlabel('x(m)')                                                      % �K�[x�b����
% %         ylabel('y(m)')                                                      % �K�[y�b����
% %                 axis([-1.5  1.5 -1.5 1.5])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         grid on                                                             % ��ܹϧΤ�����
% %    
% %         
% %         figure(3)                                                           % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(m4,m5,'r','LineWidth',1.4);                           % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('ISE of e_{att} (radian)')                                   % �K�[y�b����
% %         axis([0  50 0 1])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %              
% %         
% %         figure(4)                                                           % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(m4,m6,'r','LineWidth',1.4);                           % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('ISE of e_{pos} (m)')                                        % �K�[y�b����
% %         axis([0  50 0 0.5])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %     
% %         
% %         figure(11)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(m4,m7,'r','LineWidth',1.4);                             % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_\psi (radian)')                                           % �K�[y�b����
% %         axis([0  50 -1.1 0.5])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(12)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(m4,m8,'r','LineWidth',1.4);                           % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_\theta (radian)')                                         % �K�[y�b����
% %         axis([0  50 -0.05 0.15])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(13)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(m4,m9,'r','LineWidth',1.4);                             % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_\phi (radian)')                                           % �K�[y�b����
% %              axis([0  50 -0.02 0.1])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(14)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(m4,m10,'r','LineWidth',1.4);                               % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_z (m)')                                                   % �K�[y�b����
% %                         axis([0  50 -0.2 0.3])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(15)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(m4,m11,'r','LineWidth',1.4);                               % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_y (m)')                                                   % �K�[y�b����
% %         axis([0  50 -0.2 0.4])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% %         
% %         
% %         figure(16)                                                          % �Ыعϧε���
% %         hold on                                                             % �O�d��eø��
% %         plot(m4,m12,'r','LineWidth',1.4);                               % ø�s�ϧ�
% %         legend({'Novel IT2 N=2 ','Novel T1 N=2 ','IT2 N=2 '},'FontSize',12)                                 % �󧤼аϲK�[�Ϩ�
% %         xlabel('Time(s)')                                                   % �K�[x�b����
% %         ylabel('e_x (m)')                                                   % �K�[y�b����
% %         axis([0  50 -0.3 0.5])
% %         set(gca,'FontSize',14,'FontWeight','bold')                          % �]�w�ϧ��ݩʰѼ�
% %         box on                                                              % ��ܧ������
% % 
% %     end
% % end

%}
toc                                                                 % ��ܵ{���`����ɶ�