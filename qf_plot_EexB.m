Eex1=reshape(EexB(:,length(Dt),:),[length(Em),nE]);
nuex1=reshape(nuexB(:,length(Dt),:),[length(Em),nE]);
ndex1=reshape(ndexB(:,length(Dt),:),[length(Em),nE]);
Sz=(nuex1-ndex1)/2;
for jm=1:length(Em)
    [m, im]=sort(Sz(jm,:));
    Sz(jm,:)=m;
    Eex1(jm,:)=Eex1(jm,im);
end
cmap=colormap(linspecer(100,'sequential'));
lc=length(cmap);
for jE=1:nE
    index=fix(Sz(1,jE)/3*(lc-1))+1;
    plot(Em,Eex1(:,jE),'color',[cmap(index,:)]);
    hold on;
end
hold off;
h=colorbar;
caxis([0 3]);
set(h,'ytick',0:3);
%ylim([-5.5 0.5]);
set(gca,'FontSize',18);
set(gcf,'color','white');
set(gca,'ytick',[-9:-6]);
set(gca,'xtick',[0,0.1,0.2,0.3]);

    