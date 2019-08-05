load qf_dat_zakB_U0
phase1=abs(phase);
colormap(linspecer(100,'sequential'));
% colormap('jet');
image([min(Dt),max(Dt)],[min(Em),max(Em)],phase1,'CDataMapping','scaled');
set(gca,'Ydir','Normal','FontSize',14);
h=colorbar;
caxis([0 1]);
% xlabel('$\Delta t$','interpreter','latex');
% ylabel('$U$','interpreter','latex');
set(h,'YTick',[0:1]);
box on;
set(gcf,'color','white');
set(gca,'YTick',[0,0.1,0.2]);