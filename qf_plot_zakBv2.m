load qf_dat_zakBv2
colormap(linspecer(100,'sequential'));
% colormap('jet');
image([min(Dt),max(Dt)],[min(Em),max(Em)],phase,'CDataMapping','scaled');
set(gca,'Ydir','Normal','FontSize',14);
h=colorbar;
caxis([-0.5 0.5]);
% xlabel('$\Delta t$','interpreter','latex');
% ylabel('$U$','interpreter','latex');
set(h,'YTick',[-0.5,0,0.5]);
box on;
set(gcf,'color','white');
set(gca,'YTick',[0,0.1,0.2,0.3,0.4,0.5]);