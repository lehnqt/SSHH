Ead1=zeros(length(Ec),length(Nf));
Ead2=zeros(length(Ec),length(Nf));
pen1=zeros(length(Ec),length(Nf));
pen2=zeros(length(Ec),length(Nf));
for jn=1:length(Nf)
    for ju=1:length(Ec)
        Ead1(ju,jn)=Eadd(ju,21,jn);
        Ead2(ju,jn)=Eadd(ju,1,jn);
        pen1(ju,jn)=pen(ju,21,jn);
        pen2(ju,jn)=pen(ju,1,jn);
    end
end
Ec2=0:0.05:10;
Ead1=interp2(ones(length(Ec),1)*Nf,Ec'*ones(1,length(Nf)),Ead1,ones(length(Ec2),1)*Nf,Ec2'*ones(1,length(Nf)));
pen1=interp2(ones(length(Ec),1)*Nf,Ec'*ones(1,length(Nf)),pen1,ones(length(Ec2),1)*Nf,Ec2'*ones(1,length(Nf)));
figure;
cmap=colormap(linspecer(100,'sequential'));
lc=length(cmap);
for ju=1:length(Ec2)
   for jn=1:length(Nf)
       index=fix(pen1(ju,jn)/4*(lc-1))+1;
    scatter(Ec2(ju),Ead1(ju,jn),1,'MarkerEdgeColor',[cmap(index,:)],'MarkerFaceColor',[cmap(index,:)],'LineWidth',0.1);
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
Ead2=interp2(ones(length(Ec),1)*Nf,Ec'*ones(1,length(Nf)),Ead2,ones(length(Ec2),1)*Nf,Ec2'*ones(1,length(Nf)));
pen2=interp2(ones(length(Ec),1)*Nf,Ec'*ones(1,length(Nf)),pen2,ones(length(Ec2),1)*Nf,Ec2'*ones(1,length(Nf)));
figure;
cmap=colormap(linspecer(100,'sequential'));
lc=length(cmap);
for ju=1:length(Ec2)
   for jn=1:length(Nf)
       index=fix(pen2(ju,jn)/4*(lc-1))+1;
    scatter(Ec2(ju),Ead2(ju,jn),1,'MarkerEdgeColor',[cmap(index,:)],'MarkerFaceColor',[cmap(index,:)],'LineWidth',0.1);
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

