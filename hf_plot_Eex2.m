Eex1=reshape(Eex(:,length(Dt),:),[length(Ec),nE]);
nuex1=reshape(nuex(:,length(Dt),:),[length(Ec),nE]);
ndex1=reshape(ndex(:,length(Dt),:),[length(Ec),nE]);
Sz=(nuex1-ndex1)/2;
for jc=1:length(Ec)
    [m, im]=sort(Sz(jc,:));
    Sz(jc,:)=m;
    Eex1(jc,:)=Eex1(jc,im);
end
cmap=colormap(linspecer(100,'sequential'));
lc=length(cmap);
for jE=1:nE
    index=fix(Sz(1,jE)/6*(lc-1))+1;
    plot(Ec,Eex1(:,jE),'color',[cmap(index,:)]);
    hold on;
end
hold off;
h=colorbar;
caxis([0 6]);
set(gca,'FontSize',14);
set(gcf,'color','white');
set(h,'xtick',[0,2,4,6])
set(gca,'xtick',[0,5,10])

    