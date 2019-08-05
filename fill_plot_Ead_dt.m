Ead1=zeros(length(Dt),length(Nf));
pen1=zeros(length(Dt),length(Nf));
for jn=1:length(Nf)
    for jt=1:length(Dt)
        Ead1(jt,jn)=Eadd(41,jt,jn);
        pen1(jt,jn)=pen(41,jt,jn);
    end
end
Dt2=-0.5:0.005:0.5;
Ead1=interp2(ones(length(Dt),1)*Nf,Dt'*ones(1,length(Nf)),Ead1,ones(length(Dt2),1)*Nf,Dt2'*ones(1,length(Nf)));
pen1=interp2(ones(length(Dt),1)*Nf,Dt'*ones(1,length(Nf)),pen1,ones(length(Dt2),1)*Nf,Dt2'*ones(1,length(Nf)));
figure;
cmap=colormap(linspecer(100,'sequential'));
lc=length(cmap);
for jt=1:length(Dt2)
   for jn=1:length(Nf)
       index=fix(pen1(jt,jn)/4*(lc-1))+1;
    scatter(Dt2(jt),Ead1(jt,jn),1.5,'MarkerEdgeColor',[cmap(index,:)],'MarkerFaceColor',[cmap(index,:)],'LineWidth',0.1);
    hold on;
   end
end
h=colorbar;
caxis([0 4]);
box on;
ylim([-2.5 12]);
% set(h,'YTick',[0:4]);
set(gca,'FontSize',24);
set(gcf,'color','white');
hold off;