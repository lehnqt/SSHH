load qf_dat_zak
Ec2=Ec(2:length(Ec));
phase1=abs(phase(2:length(Ec),:));
colormap(linspecer(100,'sequential'));
% colormap('jet');
image([min(Dt),max(Dt)],[min(Ec2),max(Ec2)],phase1,'CDataMapping','scaled');
set(gca,'Ydir','Normal','FontSize',14);
h=colorbar;
caxis([0 1]);
% xlabel('$\Delta t$','interpreter','latex');
% ylabel('$U$','interpreter','latex');
set(h,'YTick',[0:1]);
box on;
set(gcf,'color','white');
%set(gca,'YTick',[0,0.1,0.2]);