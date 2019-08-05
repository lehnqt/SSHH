W0=0;
N_sites=20;
epsilon_up(1:N_sites)=0;
epsilon_down=epsilon_up;
x=1:N_sites;
Dt=-0.5:0.01:0.5;
N_up=1;
N_down=0;
erg=zeros(N_sites,length(Dt));
for jt=1:length(Dt)
    fprintf('jt = %d\n',jt);
    dt=Dt(jt);
    H0=extHbd(N_sites,N_up,N_down,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,zeros(1,N_sites),W0,x);
    eiv=eig(H0);
    erg(:,jt)=eiv;
end
for jE=1:N_sites
plot(Dt,erg(jE,:),'color',[0.9100 0.4100 0.1700],'linewidth',0.5);
hold on;
end
hold off;
% xlabel('$\Delta t/\bar{t}$','interpreter','latex');
% ylabel('$E$','interpreter','latex');
set(gca,'FontSize',14);
set(gcf,'color','white');

 dt=0.5;
 H0=extHbd(N_sites,N_up,N_down,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,zeros(1,N_sites),W0,x);
 [eif,eiv]=eig(full(H0));
 [~,nbulk]=partdens(N_sites,1,0,eif(:,N_sites/4),1);
 [site,nedge]=partdens(N_sites,1,0,eif(:,N_sites/2),1);
 
figure;
plot(site,nbulk,'d','MarkerFaceColor',[222/255 136/255 72/255],'MarkerEdgeColor',[222/255 136/255 72/255],'MarkerSize',14);
hold on
plot(site,nedge,'o','MarkerFaceColor',[68/255 114/255 196/255],'MarkerEdgeColor',[68/255 114/255 196/255],'MarkerSize',10);
set(gca,'FontSize',14);
set(gca,'xtick',[1,5,10,15,20]);
set(gca,'ytick',0:0.1:0.4);
set(gcf,'color','white');

save dat_nonint_E