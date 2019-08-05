figure;
load('dat_condctB_N12_topol.mat');
g3=g;
plot(chem/4,g3(21,:),'-r','LineWidth',1);
xlim([-20/4 60/4]);
box on;
set(gca,'FontSize',24);
set(gcf,'color','white');
hold on;
load('dat_condctB_N12_triv.mat');
g4=g;
plot(chem/4,g4(21,:),':b','LineWidth',1);
xlim([-20/4 60/4]);
box on;
set(gca,'FontSize',24);
set(gcf,'color','white');
ylim([0,2.8])
set(gca,'ydir','reverse')