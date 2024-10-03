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
    
        plot(t,xt,'r','LineWidth',1.2)
        hold all
        plot(t,ex,'b','LineWidth',1.2)

    %                plot(t,xt,'r','LineWidth',1.2)
    %set(gca,'FontWeight','bold','FontSize',16);
    hold all
    xlabel('time(sec)')
    ylabel('e_x(m)')
    legend('e_x\rm(t_k)','e_x\rm(t)')
    
    %axis([21 28-0.6 0.2])
    axis([21 28 -0.5 0.5])
    
     
     figure(3)

         plot(t,xt1,'r','LineWidth',1.2)
        hold all
         plot(t,ey,'b','LineWidth',1.2)

     %set(gca,'FontWeight','bold','FontSize',16);
     hold all
    xlabel('time(sec)')
    ylabel('e_y(m)')
    legend('e_y\rm(t_k)','e_y\rm(t)')
    axis([21 28 -0.1 0.15])

    figure(5)

    plot(t,xt2,'r','LineWidth',1.2)
    hold all
    plot(t,et,'b','LineWidth',1.2)
        %plot(t,xt2,'b','LineWidth',1.2)
    hold all
    xlabel('time(sec)')
    ylabel('e_\theta\rm(radian)')
    legend('e_\theta\rm(t_k)','e_\theta\rm(t)')
    
     axis([21 28 -0.5 0.5])

% %     figure(7)
% %     
% %         plot(t,u1,'r','LineWidth',1.2)
% %         hold all
% %     %         plot(x,y,'r','LineWidth',1.2)
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     hold all
% %     xlabel('time(sec)')
% %     ylabel('u_1')
% %     legend('u_1\rm(t_k)')
% %     axis([21 28 -1 0.5])
% %     %axis([-1.5 1.5 -0.6 1])
% %     
% %     figure(8)
% %         plot(t,u2,'r','LineWidth',1.2)
% % 
% %     %set(gca,'FontWeight','bold','FontSize',16);
% %     ylabel('u_2')
% %     xlabel('time(sec)')
% %     legend('u_2\rm(t_k)')
% %     
% %     hold all
% %     axis([21 28 -2 1.5])
end



%}


