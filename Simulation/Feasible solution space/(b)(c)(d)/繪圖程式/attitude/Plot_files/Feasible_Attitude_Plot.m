 clear
close all
clc

load Attitude_NPPF_N2
load Attitude_EPPF_N2 

figure(1)
plot(Attitude_NPPF_N2(:,1),Attitude_NPPF_N2(:,2),'ro','LineWidth',1.2,'Markersize',12)
hold on
plot(Attitude_IT2PPF_N2 (:,1),Attitude_IT2PPF_N2 (:,2),'b*','LineWidth',1,'Markersize',10)
legend('Novel SOS N=2','Existing SOS N=2','Location', 'Best')
xlabel('\gamma','fontweight','bold','fontsize',12)
ylabel('\kappa','fontweight','bold','fontsize',12)
axis([0.1 1 1 10])
set(gca,'fontweight','bold','fontsize',12)
box on
hold off