colormap(linspecer(100,'sequential'));
EmL2=0:0.001:1;
dat=interp2(chem.*ones(length(EmL),1),EmL'.*ones(1,length(chem)),g,chem.*ones(length(EmL2),1),EmL2'.*ones(1,length(chem)));
% colormap('jet');
image([chem(1)/4,chem(length(chem))/4],[EmL(1)/4,EmL(length(EmL))/4],dat,'CDataMapping','scaled');

h=colorbar;
% xlabel('$\Delta t$','interpreter','latex');
% ylabel('$U$','interpreter','latex');
box on;
set(gca,'Ydir','Normal','FontSize',24);
set(gcf,'color','white');
xlim([-5/4,5/4]);
ylim([0,1/4]);
set(h,'FontSize',24)