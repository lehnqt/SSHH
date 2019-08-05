Eex1=reshape(Eex(length(Ec),:,:),[length(Dt),nE]);
nuex1=reshape(nuex(length(Ec),:,:),[length(Dt),nE]);
ndex1=reshape(ndex(length(Ec),:,:),[length(Dt),nE]);
Sz=(nuex1-ndex1)/2;
for jt=1:length(Dt)
    [m, im]=sort(Sz(jt,:));
    Sz(jt,:)=m;
    Eex1(jt,:)=Eex1(jt,im);
end
cmap=colormap(linspecer(100,'sequential'));
lc=length(cmap);
for jE=1:nE
    index=fix(Sz(1,jE)/3*(lc-1))+1;
    plot(Dt,Eex1(:,jE),'color',[cmap(index,:)]);
    hold on;
end
hold off;
h=colorbar;
caxis([0 3]);
set(h,'ytick',0:3);
%ylim([-5.5 0.5]);
set(gca,'FontSize',14);
set(gcf,'color','white');

    