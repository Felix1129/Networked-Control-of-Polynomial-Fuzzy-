 clear 
clc
close all
jj=1
for jj=1
data = csvread('TH4N_N2.csv');
% t = rmmissing(data(:,1));
xr = rmmissing(data(:,1));
yr = rmmissing(data(:,2));
thr = rmmissing(data(:,3));
xc = rmmissing(data(:,4));
yc = rmmissing(data(:,5));
thc = rmmissing(data(:,6));
ex = rmmissing(data(:,7));
ey = rmmissing(data(:,8));
eth = rmmissing(data(:,9));
 en = rmmissing(data(:,11));
th=1;
t = 1:length(ex);
figure(1)
if th==1
  plot(xr,yr,'k--',xc,yc,'r')
end
if th==2
 plot(xc,yc,'b')
end
if th==3
 plot(xc,yc,'k')
end
legend('Reference','Novel SOS N=2','Existing SOS N=2')

axis([-1 4.5  -1 4.5])
xlabel('{\itx}(m)')
ylabel('{\ity}(m)')
hold on

figure(3)
if th==1
plot(t/55,ex,'r')
end
if th==2
plot(t/55,ex,'b')
end
if th==3
plot(t/55,ex,'k')
end
xlabel('time(sec)')
axis([0 40 -1.5 1.5])
ylabel('{\ite_x}(m)')
legend('Novel SOS N=2','Existing SOS N=2')

hold on

figure(5)
if th==1
plot(t/55,ey,'r')
end
if th==2
plot(t/55,ey,'b')
end
if th==3
plot(t/55,ey,'k')
end
axis([0 40 -1 1])
xlabel('time(sec)')
ylabel('{\ite_y}(m)')
legend('Novel SOS N=2','Existing SOS N=2')

hold on

figure(7)
if th==1
plot(t/55,eth,'r')
end
if th==2
plot(t/55,eth,'b')
end
if th==3
plot(t/55,eth,'k')
end
axis([0 40 -1 1])
xlabel('time(sec)')
ylabel('{\ite_\theta}(radian)')
legend('Novel SOS N=2','Existing SOS N=2')

hold on

figure(8)
if th==1
plot(t/55,en,'r')
end
if th==2
plot(t/55,en,'b')
end
if th==3
plot(t/55,en,'k')
end
axis([0 40 0 2])
xlabel('time(sec)')
ylabel('{\iten}({\itt})')
legend('Novel SOS N=2','Existing SOS N=2')

hold on


end

for jj=1
data = csvread('TH4_N2.csv');
% t = rmmissing(data(:,1));
xr = rmmissing(data(:,1));
yr = rmmissing(data(:,2));
thr = rmmissing(data(:,3));
xc = rmmissing(data(:,4));
yc = rmmissing(data(:,5));
thc = rmmissing(data(:,6));
ex = rmmissing(data(:,7));
ey = rmmissing(data(:,8));
eth = rmmissing(data(:,9));
 en = rmmissing(data(:,11));
th=2;
t = 1:length(ex);
figure(1)
if th==1
  plot(xr,yr,'k--',xc,yc,'r')
end
if th==2
 plot(xc,yc,'b')
end
if th==3
 plot(xc,yc,'k')
end
legend('Reference','Novel SOS N=2','Existing SOS N=2')

axis([-1 4.5  -1 4.5])
xlabel('{\itx}(m)')
ylabel('{\ity}(m)')
hold on

figure(3)
if th==1
plot(t/55,ex,'r')
end
if th==2
plot(t/55,ex,'b')
end
if th==3
plot(t/55,ex,'k')
end
xlabel('time(sec)')
axis([0 40 -1.5 1.5])
ylabel('{\ite_x}(m)')
legend('Novel SOS N=2','Existing SOS N=2')

hold on

figure(5)
if th==1
plot(t/55,ey,'r')
end
if th==2
plot(t/55,ey,'b')
end
if th==3
plot(t/55,ey,'k')
end
axis([0 40 -1 1])
xlabel('time(sec)')
ylabel('{\ite_y}(m)')
legend('Novel SOS N=2','Existing SOS N=2')

hold on

figure(7)
if th==1
plot(t/55,eth,'r')
end
if th==2
plot(t/55,eth,'b')
end
if th==3
plot(t/55,eth,'k')
end
axis([0 40 -1 1])
xlabel('time(sec)')
ylabel('{\ite_\theta}(radian)')
legend('Novel SOS N=2','Existing SOS N=2')

hold on

figure(8)
if th==1
plot(t/55,en,'r')
end
if th==2
plot(t/55,en,'b')
end
if th==3
plot(t/55,en,'k')
end
axis([0 40 0 2])
xlabel('time(sec)')
ylabel('{\iten}({\itt})')
legend('Novel SOS N=2','Existing SOS N=2')

hold on


end


