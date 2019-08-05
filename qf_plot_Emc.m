load qf_dat_EexB_Ec_Dt_kara.mat
nu=nuexB(:,:,:,1);
Emc=zeros(length(Dt),length(Ec));
for jt=1:length(Dt)
    for jc=1:length(Ec)
    v=nu(:,jt,jc);
    ind=find(v==N_fill,1,'first');
    Emc(jt,jc)=Em(ind);
    end
end
colormap(linspecer(100,'sequential'));
% colormap('jet');
dat=interp2(Emc',4);
image([min(Dt),max(Dt)],[min(Ec),max(Ec)],dat,'CDataMapping','scaled');
caxis([0 1]);
set(gca,'Ydir','Normal','FontSize',24);
h=colorbar;
% xlabel('$\Delta t$','interpreter','latex');
% ylabel('$U$','interpreter','latex');
set(h,'YTick',[0.1,0.5,0.9]);
box on;
set(gcf,'color','white');
    