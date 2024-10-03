 close all
clc

load Altitude_NPPF_N2
load Altitude_EPPF_N2

figure(1)
plot(Altitude_NPPF_N2(:,1),Altitude_NPPF_N2(:,2),'ro','LineWidth',1.2,'Markersize',12)
hold on
plot(Altitude_EPPF_N2(:,1),Altitude_EPPF_N2(:,2),'b*','LineWidth',1,'Markersize',10)
legend('Novel SOS N=2','Existing SOS N=2','Location', 'Best')
xlabel('\gamma','fontweight','bold','fontsize',12)
ylabel('\kappa','fontweight','bold','fontsize',12)
axis([0.1 1 1 10])
set(gca,'fontweight','bold','fontsize',12)
box on
hold off