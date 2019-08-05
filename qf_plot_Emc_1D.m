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
plot(Dt,Emc(:,1),'-b','LineWidth',1);
set(gca,'fontsize',24);
ylim([0,1]);
set(gcf,'color','white');
figure;
plot(Ec,Emc(length(Dt),:),'-b','LineWidth',1);
set(gca,'fontsize',24,'ytick',[]);
ylim([0,1]);
set(gcf,'color','white');