 clear
close all
clc








for gg=1

N1 = xlsread('N_N2.xlsx');
xr_1 = rmmissing(N1(:,1));
yr_1 = rmmissing(N1(:,5));
zr_1 = rmmissing(N1(:,9));

xc_1 = rmmissing(N1(:,2));
yc_1 = rmmissing(N1(:,6));
zc_1 = rmmissing(N1(:,10));

deg1r_1=rmmissing(N1(:,13));
deg2r_1=rmmissing(N1(:,17));
deg3r_1=rmmissing(N1(:,21));

deg1c_1=rmmissing(N1(:,14));
deg2c_1=rmmissing(N1(:,18));
deg3c_1=rmmissing(N1(:,22));

t_1 = 1:length(xr_1);

for i = 1:t_1(end)
    ex_1(i) = xc_1(i)-xr_1(i);
    ey_1(i) = yc_1(i)-yr_1(i);
    ez_1(i) = zc_1(i)-zr_1(i);
end
for i = 1:t_1(end)
    edeg1_1(i) = (deg1c_1(i)-deg1r_1(i))*3.14159/180;
    edeg2_1(i) = (deg2c_1(i)-deg2r_1(i))*3.14159/180;
    edeg3_1(i) = (deg3c_1(i)-deg3r_1(i))*3.14159/180;
end
for i = 1:t_1(end)
   en_a(i) =  (ex_1(i))^2 + (ey_1(i))^2 + (ez_1(i))^2;
end
for i = 1:t_1(end)
   en_b(i) =  (edeg1_1(i))^2 + (edeg2_1(i))^2 + (edeg3_1(i))^2;
end
en_pos_1 = cumsum(en_a);
en_att_1 = cumsum(en_b);
end
















%бе╕Б
for gg=1
N2 = xlsread('N2.xlsx');
xr_2 = rmmissing(N2(:,1));
yr_2 = rmmissing(N2(:,5));
zr_2 = rmmissing(N2(:,9));

xc_2 = rmmissing(N2(:,2));
yc_2 = rmmissing(N2(:,6));
zc_2 = rmmissing(N2(:,10));

deg1r_2=rmmissing(N2(:,13));
deg2r_2=rmmissing(N2(:,17));
deg3r_2=rmmissing(N2(:,21));

deg1c_2=rmmissing(N2(:,14));
deg2c_2=rmmissing(N2(:,18));
deg3c_2=rmmissing(N2(:,22));

t_2 = 1:length(xr_2);

for i = 1:t_2(end)
    ex_2(i) = xc_2(i)-xr_2(i);
    ey_2(i) = yc_2(i)-yr_2(i);
    ez_2(i) = zc_2(i)-zr_2(i);
end
for i = 1:t_2(end)
    edeg1_2(i) = (deg1c_2(i)-deg1r_2(i))*3.14159/180;
    edeg2_2(i) = (deg2c_2(i)-deg2r_2(i))*3.14159/180;
    edeg3_2(i) = (deg3c_2(i)-deg3r_2(i))*3.14159/180;
end

for i = 1:t_2(end)
   en_c(i) =  (ex_2(i))^2 + (ey_2(i))^2 + (ez_2(i))^2;
end
for i = 1:t_2(end)
   en_d(i) =  (edeg1_2(i))^2 + (edeg2_2(i))^2 + (edeg3_2(i))^2;
end
en_pos_2 = cumsum(en_c);
en_att_2 = cumsum(en_d);
end




%}




%}

%%

figure(1)
hold on
plot(t_1,ex_1,'r','LineWidth',1.4)
plot(t_2,ex_2,'b','LineWidth',1.4)


xlabel('Time(s)')
ylabel('e_x(m)')
axis([0 720 -3 4])
xticks([0 180 360 540 720 900])
xticklabels({'0','20','40','60','80','100'})
legend({'Novel SOS N=2','Existing SOS N=2'},'FontSize',12)
set(gca,'FontSize',14,'FontWeight','bold')
box on

figure(2)
hold on
plot(t_1,ey_1,'r','LineWidth',1.4)
plot(t_2,ey_2,'b','LineWidth',1.4)

xlabel('Time(s)')
ylabel('e_y(m)')
axis([0 720 -3 3])
xticks([0 150 400 600 800 1000])
xticklabels({'0','20','40','60','80','100'})
legend({'Novel SOS N=2','Existing SOS N=2'},'FontSize',12)
set(gca,'FontSize',14,'FontWeight','bold')
box on

figure(3)
hold on
plot(t_1,ez_1,'r','LineWidth',1.4)
plot(t_2,ez_2,'b','LineWidth',1.4)
xlabel('Time(s)')
ylabel('e_z(m)')
axis([0 720 -1 2])
xticks([0 180 360 540 720 900])
xticklabels({'0','20','40','60','80','100'})
legend({'Novel SOS N=2','Existing SOS N=2'},'FontSize',12)
set(gca,'FontSize',14,'FontWeight','bold')
box on

figure(4)
hold on

plot(t_1,en_pos_1/9,'r','LineWidth',1.4)
plot(t_2,en_pos_2/9,'b','LineWidth',1.4)

xlabel('Time(s)')
ylabel('ISE_{pos}')
axis([0 720 0 300])
xticks([0 180 360 540 720 900])
xticklabels({'0','20','40','60','80','100'})
legend({'Novel SOS N=2','Existing SOS N=2'},'FontSize',12)
set(gca,'FontSize',14,'FontWeight','bold')
box on

figure(5)
hold on
plot(xr_1(1:end-250),yr_1(1:end-250),'k--','LineWidth',1.4)
plot(xc_1(1:end-250),yc_1(1:end-250),'r','LineWidth',1.4)
plot(xc_2(1:end-250),yc_2(1:end-250),'b','LineWidth',1.4)


xlabel('x(m)')
ylabel('y(m)')
 axis([-12 12 -7 12 ])
legend({'Reference','Novel SOS N=2','Existing SOS N=2'},'FontSize',10)
set(gca,'FontSize',14,'FontWeight','bold')
box on

figure(6)
hold on
plot3(xr_1(1:end-250),yr_1(1:end-250),zr_1(1:end-250),'k--','LineWidth',1.4)
plot3(xc_1(1:end-250),yc_1(1:end-250),zc_1(1:end-250),'r','LineWidth',1.4)
plot3(xc_2(1:end-250),yc_2(1:end-250),zc_2(1:end-250),'b','LineWidth',1.4)


grid on
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
view(-15,20);
% axis([-11 2 -6 6 0 7])
legend({'Reference','Novel IT2 N=2','Novel T1 N=2','IT2 N=2'},'FontSize',10)
set(gca,'FontSize',14,'FontWeight','bold')
box on


figure(7)
hold on
plot(t_1,edeg1_1,'r','LineWidth',1.4)
plot(t_2,edeg1_2,'b','LineWidth',1.4)


xlabel('Time(s)')
ylabel('e_\phi(m)')
axis([0 720 -0.5 0.5])
xticks([0 180 360 540 720 900])
xticklabels({'0','20','40','60','80','100'})
legend({'Novel SOS N=2','Existing SOS N=2'},'FontSize',12)
set(gca,'FontSize',14,'FontWeight','bold')
box on


figure(8)
hold on
plot(t_1,edeg2_1,'r','LineWidth',1.4)
plot(t_2,edeg2_2,'b','LineWidth',1.4)


xlabel('Time(s)')
ylabel('e_\theta(m)')
axis([0 720 -0.5 0.5])
xticks([0 180 360 540 720 900])
xticklabels({'0','20','40','60','80','100'})
legend({'Novel SOS N=2','Existing SOS N=2'},'FontSize',12)
set(gca,'FontSize',14,'FontWeight','bold')
box on


figure(9)
hold on
plot(t_1,edeg3_1,'r','LineWidth',1.4)
plot(t_2,edeg3_2,'b','LineWidth',1.4)


xlabel('Time(s)')
ylabel('e_\psi(m)')
axis([0 720 -1 1])
xticks([0 180 360 540 720 900])
xticklabels({'0','20','40','60','80','100'})
legend({'Novel SOS N=2','Existing SOS N=2'},'FontSize',12)
set(gca,'FontSize',14,'FontWeight','bold')
box on



figure(10)
hold on

plot(t_1,en_att_1/9,'r','LineWidth',1.4)
plot(t_2,en_att_2/9,'b','LineWidth',1.4)

xlabel('Time(s)')
ylabel('ISE_{att}')
axis([0 720 0 2])
xticks([0 180 360 540 720 900])
xticklabels({'0','20','40','60','80','100'})
legend({'Novel SOS N=2','Existing SOS N=2'},'FontSize',12)
set(gca,'FontSize',14,'FontWeight','bold')
box on





