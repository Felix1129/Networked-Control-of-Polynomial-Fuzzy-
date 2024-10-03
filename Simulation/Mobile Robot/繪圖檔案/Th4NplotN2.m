 clear all
close all
clc
tic
double exmax
global max_ex_in_ex max_ey_in_ey max_ey_and_ex
global col_ex row_ex col_ey row_ey
global m_in_i1 m_in_i2 m_in_i3 m_in_i4 m_in_i5 m_in_i6
global g flag_hold_dis
global i1 i2 i3 i4 i5 i6
global reg_sum1 en1 en2 en3
global reg1 reg2 reg3 reg4 reg5 reg6 reg7
global xdis xdis2 xdis3 ydis ydis2 ydis3
global flag
eyfin_max=0;
eymax = -1000;
eyreg = 1000;
exfin_max=0;
exmax = -1000;
exreg = 1000;

coun=0;
ex_dif=-1000;
sum_dif=0;
last_sum_dif=100000000;

flag_hold_dis=0;
the=1;
reg_sum1=0;
%
for gg=1

                            flag=1;
                            sim('NSOSN2.mdl');  %%家览郎
                            [col_ex,row_ex]=size(ex);
                            [col_ey,row_ey]=size(ey);
                            [col_x,row_x]=size(x);
                            [col_y,row_y]=size(y);
                            for i =1: col_ex
                                  if flag_dis(i,1) ~= 0 && flag_hold_dis==0
                                    flag_hold_dis=1;
                                    xdis= x(i,1);
                                    ydis= y(i,1);
                                end
                                reg_sum1=(ex(i,1)*ex(i,1)+ey(i,1)*ey(i,1))/4077+reg_sum1;
                                en1(i,1)=reg_sum1;
                                if exmax < ex(i,1)*1000 ;
                                    exmax = ex(i,1)*1000;
                                    max_ex_in_ex =i;
                                end%%程禬禫exmax
                                if eymax < ey(i,1)*1000 ;
                                    eymax = ey(i,1)*1000;
                                    max_ey_in_ey =i;
                                end%%程禬禫eymax
                                %xr_new-x,yr_new-yr璶―ゑ顶琿ぶ
                                y_dif(i,1)= ey(i,1)-y(i,1);
                                x_dif(i,1)= ex(i,1)-x(i,1);
                                all_dif(i,1)=sqrt(y_dif(i,1)*y_dif(i,1)+x_dif(i,1)*x_dif(i,1));
                                sum_dif=sum_dif+all_dif(i,1);
                            end
                            exfin_max = exmax;
                            eyfin_max = eymax;
                            exmax=-1000;%%程禬禫耴-1000
                            eymax=-1000;%%程禬禫耴-1000
                            sum_dif=0;
    
    toc
  %舼︹絬瓜程よ–Ωメアぃ
  %璶癘材Ω禲戈
    if the == 1
        reg1=t;
        reg2=ex;
        reg3=ey;
        reg4=et;
        reg5=x;
        reg6=y;
        reg7=en1;       
    end
    
    figure(1)
    if the == 1
        plot(t,ex,'r','LineWidth',1.2)
    end
    if the == 2
        plot(t,ex,'b','LineWidth',1.2)
    end
%     if the == 3
%         plot(t,ex,'g','LineWidth',1.2)
%     end
    %                plot(t,xt,'r','LineWidth',1.2)
    %set(gca,'FontWeight','bold','FontSize',16);
    hold all
    xlabel('time (sec)')
    ylabel('e_x')
    legend('Novel SOS N=2','Existing SOS N=2')
    
    %axis([0 40-0.6 0.2])
    axis([0 40 -1.2 1])
    
    
    figure(3)
    if the ==1
        plot(t,ey,'r','LineWidth',1.2)
    end
    if the ==2
        plot(t,ey,'b','LineWidth',1.2)
    end
%     if the ==3
%         plot(t,ey,'g','LineWidth',1.2)
%     end
    %set(gca,'FontWeight','bold','FontSize',16);
    hold all
    xlabel('time (sec)')
    ylabel('e_y')
    legend('Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
    axis([0 40 -1 1])
    %axis([0 40-0.5 0.3])
    
    figure(5)
    if the ==1
        plot(t,et,'r','LineWidth',1.2)
    end
    if the ==2
        plot(t,et,'b','LineWidth',1.2)
        %plot(t,xt2,'b','LineWidth',1.2)
    end
%     if the ==3
%         plot(t,et,'g','LineWidth',1.2)
%     end
    %set(gca,'FontWeight','bold','FontSize',16);
    hold all
    xlabel('time (sec)')
    ylabel('e_\theta\rm(radian)')
    legend('Novel SOS N=2','Existing SOS N=2')
    
    axis([0 40 -2 2])
    %礚axis
    figure(7)
    

    if the ==1
        plot(xr,yr,'k--',x,y,'r','LineWidth',1.2)
        
    end
    if the ==2
        plot(x,y,'b','LineWidth',1.2)
    end
%     if the ==3
%         plot(x,y,'g','LineWidth',1.2)
%     end
    %         plot(x,y,'r','LineWidth',1.2)
    %set(gca,'FontWeight','bold','FontSize',16);
    hold all
    xlabel('x\rm (m)')
    ylabel('y\rm (m)')
    legend('Reference','Novel SOS N=2','Existing SOS N=2')
    axis([ -1 4 -1 4])
    %axis([-1.5 1.5 -0.6 1])
    
    figure(8)
    if the ==1
        plot(t,en1,'r','LineWidth',1.2)
    end
    if the ==2
        plot(t,en1,'b','LineWidth',1.2)
    end
    if the ==3
        plot(t,en1,'g','LineWidth',1.2)
    end
    %set(gca,'FontWeight','bold','FontSize',16);
    ylabel('ISE')
    xlabel('time(sec)')
    legend('Novel SOS N=2','Existing SOS N=2')
    
    hold all
    axis([0 40 0 2])
    %    axis([0 80 0 2])
end

% % the=2;
% % reg_sum1=0;
% % flag_hold_dis=0;
% % 
% % for gg=1
% % 
% %                             flag=2;
% %                             sim('NIT2SOSN2.mdl');  %%家览郎
% % 
% %                             [col_ex,row_ex]=size(ex);
% %                             [col_ey,row_ey]=size(ey);
% %                             [col_x,row_x]=size(x);
% %                             [col_y,row_y]=size(y);
% %                             for i =1: col_ex
% %                                 if flag_dis(i,1) ~= 0 && flag_hold_dis==0
% %                                     flag_hold_dis=1;
% %                                     xdis2= x(i,1);
% %                                     ydis2= y(i,1);
% %                                 end
% %                                 reg_sum1=(ex(i,1)*ex(i,1)+ey(i,1)*ey(i,1))/4077+reg_sum1;
% %                                 en2(i,1)=reg_sum1;
% %                                 if exmax < ex(i,1)*1000 ;
% %                                     exmax = ex(i,1)*1000;
% %                                     max_ex_in_ex =i;
% %                                 end%%程禬禫exmax
% %                                 if eymax < ey(i,1)*1000 ;
% %                                     eymax = ey(i,1)*1000;
% %                                     max_ey_in_ey =i;
% %                                 end%%程禬禫eymax
% %                                 %xr_new-x,yr_new-yr璶―ゑ顶琿ぶ
% %                                 y_dif(i,1)= ey(i,1)-y(i,1);
% %                                 x_dif(i,1)= ex(i,1)-x(i,1);
% %                                 all_dif(i,1)=sqrt(y_dif(i,1)*y_dif(i,1)+x_dif(i,1)*x_dif(i,1));
% %                                 sum_dif=sum_dif+all_dif(i,1);
% %                             end
% %                          
% %                             exfin_max = exmax;
% %                             eyfin_max = eymax;
% %                             exmax=-1000;%%程禬禫耴-1000
% %                             eymax=-1000;%%程禬禫耴-1000
% %                             sum_dif=0;
% %  
% %     
% %     toc
% %   
% %     
% %     
% %     figure(1)
% %     if the == 1
% %         plot(t,ex,'r','LineWidth',1.2)
% %     end
% %     if the == 2
% %         plot(t,ex,'b','LineWidth',1.2)
% %     end
% %     if the == 3
% %         plot(t,ex,'g','LineWidth',1.2)
% %     end
% %     %                plot(t,xt,'r','LineWidth',1.2)
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     hold all
% %     xlabel('time(sec)')
% %     ylabel('e_x')
% %     legend('Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     
% %     %axis([0 40-0.6 0.2])
% %     axis([0 40 -1.2 1])
% %     
% %     
% %     figure(3)
% %     if the ==1
% %         plot(t,ey,'r','LineWidth',1.2)
% %     end
% %     if the ==2
% %         plot(t,ey,'b','LineWidth',1.2)
% %     end
% %     if the ==3
% %         plot(t,ey,'g','LineWidth',1.2)
% %     end
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     hold all
% %     xlabel('time(sec)')
% %     ylabel('e_y')
% %     legend('Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     axis([0 40 -1 1])
% %     %axis([0 40-0.5 0.3])
% %     
% %     figure(5)
% %     if the ==1
% %         plot(t,et,'r','LineWidth',1.2)
% %     end
% %     if the ==2
% %         plot(t,et,'b','LineWidth',1.2)
% %         %plot(t,xt2,'b','LineWidth',1.2)
% %     end
% %     if the ==3
% %         plot(t,et,'g','LineWidth',1.2)
% %     end
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     hold all
% %     xlabel('time(sec)')
% %     ylabel('e_\theta\rm(rad)')
% %     legend('Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     
% %     axis([0 40 -2 2])
% %     %礚axis
% %     figure(7)
% %     
% %     %
% %     % if the ==1
% %     %     plot(xr,yr,'k--','LineWidth',1.2)
% %     %     hold on
% %     %     plot(x,y,'r','LineWidth',1.2)
% %     % end
% %     % %
% %     if the ==1
% %         plot(xr,yr,'k--',x,y,'r','LineWidth',1.2)
% %         
% %     end
% %     if the ==2
% %         plot(x,y,'b','LineWidth',1.2)
% %     end
% %     if the ==3
% %         plot(x,y,'g','LineWidth',1.2)
% %     end
% %     %         plot(x,y,'r','LineWidth',1.2)
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     hold all
% %     xlabel('x\rm(m)')
% %     ylabel('y\rm(m)')
% %     legend('Reference','Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     axis([ -1 4 -1 4])
% %     %axis([-1.5 1.5 -0.6 1])
% %     
% %     figure(8)
% %     if the ==1
% %         plot(t,en2,'r','LineWidth',1.2)
% %     end
% %     if the ==2
% %         plot(t,en2,'b','LineWidth',1.2)
% %     end
% %     if the ==3
% %         plot(t,en2,'g','LineWidth',1.2)
% %     end
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     ylabel('ISE\rm(t\rm)')
% %     xlabel('time(sec)')
% %     legend('Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     
% %     hold all
% %     axis([0 40 0 2])
% %     %    axis([0 80 0 2])
% % end
% % 
% % the=3;
% % reg_sum1=0;
% % %
% % for gg=1
% %                             flag=3;
% %                             sim('NIT2SOSN2.mdl');  %%家览郎
% %                             [col_ex,row_ex]=size(ex);
% %                             [col_ey,row_ey]=size(ey);
% %                             [col_x,row_x]=size(x);
% %                             [col_y,row_y]=size(y);
% %                             for i =1: col_ex
% %                                 if flag_dis(i,1) ~= 0 && flag_hold_dis==0
% %                                     flag_hold_dis=1;
% %                                     xdis3= x(i,1);
% %                                     ydis3= y(i,1);
% %                                 end
% %                                 reg_sum1=(ex(i,1)*ex(i,1)+ey(i,1)*ey(i,1))/4077+reg_sum1;
% %                                 en3(i,1)=reg_sum1;
% %                                 if exmax < ex(i,1)*1000 ;
% %                                     exmax = ex(i,1)*1000;
% %                                     max_ex_in_ex =i;
% %                                 end%%程禬禫exmax
% %                                 if eymax < ey(i,1)*1000 ;
% %                                     eymax = ey(i,1)*1000;
% %                                     max_ey_in_ey =i;
% %                                 end%%程禬禫eymax
% %                                 %xr_new-x,yr_new-yr璶―ゑ顶琿ぶ
% %                                 y_dif(i,1)= ey(i,1)-y(i,1);
% %                                 x_dif(i,1)= ex(i,1)-x(i,1);
% %                                 all_dif(i,1)=sqrt(y_dif(i,1)*y_dif(i,1)+x_dif(i,1)*x_dif(i,1));
% %                                 sum_dif=sum_dif+all_dif(i,1);
% %                             end
% %                             
% %   
% % 
% %     toc
% %     %end
% %     
% %     figure(1)
% %     if the == 1
% %         plot(t,ex,'r','LineWidth',1.2)
% %     end
% %     if the == 2
% %         plot(t,ex,'b','LineWidth',1.2)
% %     end
% %     if the == 3
% %         plot(t,ex,'g','LineWidth',1.2)
% %     end
% %     %                plot(t,xt,'r','LineWidth',1.2)
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     hold all
% %     xlabel('time(sec)')
% %     ylabel('e_x')
% %     legend('Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     
% %     %axis([0 40-0.6 0.2])
% %     axis([0 40 -1.2 1])
% %     
% %     
% %     figure(3)
% %     if the ==1
% %         plot(t,ey,'r','LineWidth',1.2)
% %     end
% %     if the ==2
% %         plot(t,ey,'b','LineWidth',1.2)
% %     end
% %     if the ==3
% %         plot(t,ey,'g','LineWidth',1.2)
% %     end
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     hold all
% %     xlabel('time(sec)')
% %     ylabel('e_y')
% %     legend('Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     axis([0 40 -1 1])
% %     %axis([0 40-0.5 0.3])
% %     
% %     figure(5)
% %     if the ==1
% %         plot(t,et,'r','LineWidth',1.2)
% %     end
% %     if the ==2
% %         plot(t,et,'b','LineWidth',1.2)
% %         %plot(t,xt2,'b','LineWidth',1.2)
% %     end
% %     if the ==3
% %         plot(t,et,'g','LineWidth',1.2)
% %     end
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     hold all
% %     xlabel('time(sec)')
% %     ylabel('e_\theta\rm(rad)')
% %     legend('Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     
% %     axis([0 40 -2 2])
% %     %礚axis
% %     figure(7)
% %     
% %     %
% %     % if the ==1
% %     %     plot(xr,yr,'k--','LineWidth',1.2)
% %     %     hold on
% %     %     plot(x,y,'r','LineWidth',1.2)
% %     % end
% %     % %
% %     if the ==1
% %         plot(xr,yr,'k--',x,y,'r','LineWidth',1.2)
% %         
% %     end
% %     if the ==2
% %         plot(x,y,'b','LineWidth',1.2)
% %     end
% %     if the ==3
% %         plot(x,y,'g','LineWidth',1.2)
% %     end
% %     %         plot(x,y,'r','LineWidth',1.2)
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     hold all
% %     xlabel('x\rm(m)')
% %     ylabel('y\rm(m)')
% %     legend('Reference','Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     axis([ -1 4 -1 4])
% %     %axis([-1.5 1.5 -0.6 1])
% %     
% %     figure(8)
% %     if the ==1
% %         plot(t,en3,'r','LineWidth',1.2)
% %     end
% %     if the ==2
% %         plot(t,en3,'b','LineWidth',1.2)
% %     end
% %     if the ==3
% %         plot(t,en3,'g','LineWidth',1.2)
% %     end
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     ylabel('ISE\rm(t\rm)')
% %     xlabel('time(sec)')
% %     legend('Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     
% %     hold all
% %     axis([0 40 0 2])
% %     %    axis([0 80 0 2])
% % end
% % 
% % 
% % 
% % 
% % 
% % the=1
% %  
% %     figure(1)
% %     if the == 1
% %         plot(reg1,reg2,'r','LineWidth',1.2)
% %     end
% %     if the == 2
% %         plot(t,ex,'b','LineWidth',1.2)
% %     end
% %     if the == 3
% %         plot(t,ex,'g','LineWidth',1.2)
% %     end
% %     %                plot(t,xt,'r','LineWidth',1.2)
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     hold all
% %     xlabel('time(sec)')
% %     ylabel('e_x\rm(m)')
% %     legend('Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     
% %     %axis([0 40-0.6 0.2])
% %     axis([0 40 -0.6 0.4])
% %     
% %     
% %     figure(3)
% %     if the ==1
% %         plot(reg1,reg3,'r','LineWidth',1.2)
% %     end
% %     if the ==2
% %         plot(t,ey,'b','LineWidth',1.2)
% %     end
% %     if the ==3
% %         plot(t,ey,'g','LineWidth',1.2)
% %     end
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     hold all
% %     xlabel('time(sec)')
% %     ylabel('e_y\rm(m)')
% %     legend('Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     axis([0 40 -0.7 0.4])
% %     %axis([0 40-0.5 0.3])
% %     
% %     figure(5)
% %     if the ==1
% %         plot(reg1,reg4,'r','LineWidth',1.2)
% %     end
% %     if the ==2
% %         plot(t,et,'b','LineWidth',1.2)
% %         %plot(t,xt2,'b','LineWidth',1.2)
% %     end
% %     if the ==3
% %         plot(t,et,'g','LineWidth',1.2)
% %     end
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     hold all
% %     xlabel('time(sec)')
% %     ylabel('e_\theta\rm(radian)')
% %     legend('Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     
% %     axis([0 40 -0.8 1.2])
% %     %礚axis
% %     figure(7)
% %     
% %     %
% %     % if the ==1
% %     %     plot(xr,yr,'k--','LineWidth',1.2)
% %     %     hold on
% %     %     plot(x,y,'r','LineWidth',1.2)
% %     % end
% %     % %
% %     if the ==1
% %         plot(reg5,reg6,'r','LineWidth',1.2)
% %         
% %     end
% %     if the ==2
% %         plot(reg5,reg6,'b','LineWidth',1.2)
% %     end
% %     if the ==3
% %         plot(reg5,reg6,'g','LineWidth',1.2)
% %     end
% %     %         plot(x,y,'r','LineWidth',1.2)
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     hold all
% %     
% %     aa=0;
% %     bb=0;
% %     plot(aa,bb,'ro','MarkerFaceColor','r')
% %      text_init={'initial\rightarrow ';'point'}
% %     text(0,0,text_init,'HorizontalAlignment','left')
% %     xlabel('x\rm(m)')
% %     ylabel('y\rm(m)')
% %    
% %     plot(xdis3,ydis3,'g.','LineWidth',1.2,'Markersize',15)
% %     hold on
% %     
% %     plot(xdis2,ydis2,'b.','LineWidth',1.2,'Markersize',15)
% %     hold on
% %     
% %     plot(xdis,ydis,'r.','LineWidth',1.2,'Markersize',15)
% %     hold on
% %             text_h={'disturbance\rightarrow ';' happened'}
% %         text(xdis ,ydis, text_h);
% %     axis([ -1 3 -1 3])
% %     legend('Reference','Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     %axis([-1.5 1.5 -0.6 1])
% %     
% %     figure(8)
% %     if the ==1
% %         plot(reg1,reg7,'r','LineWidth',1.2)
% %     end
% %     if the ==2
% %         plot(reg1,reg7,'b','LineWidth',1.2)
% %     end
% %     if the ==3
% %         plot(reg1,reg7,'g','LineWidth',1.2)
% %     end
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     ylabel('ISE')
% %     xlabel('time(sec)')
% %     legend('Novel IT2 N=2','Novel T1 N=2','IT2 N=2')
% %     
% %     hold all
% %     axis([0 40 0 1.5])
    %    axis([0 80 0 2])

%}


