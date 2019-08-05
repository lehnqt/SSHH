pen1=zeros(length(Ec),length(Nf));
pen2=zeros(length(Ec),length(Nf));
for jn=1:length(Nf)
    for ju=1:length(Ec)
        pen1(ju,jn)=pen(ju,21,jn);
        pen2(ju,jn)=pen(ju,1,jn);
    end
end
Nf1=ones(length(Ec),1)*Nf;
Ec1=Ec'*ones(1,length(Nf));
Ec2=0:0.01:10;
Nf2=ones(length(Ec2),1)*Nf;
Ec2=Ec2'*ones(1,length(Nf));
pen1=interp2(Nf1,Ec1,pen1,Nf2,Ec2);
pen2=interp2(Nf1,Ec1,pen2,Nf2,Ec2);
fig1=figure;
mesh(Nf2,Ec2,pen1);
xlim([0 25]);
colormap(linspecer(100,'sequential'));
box on;
xticks([2:2:24]);
set(gca,'FontSize',12);
set(gcf,'color','white');
fig2=figure;
mesh(Nf2,Ec2,pen2);
xlim([0 25]);
colormap(linspecer(100,'sequential'));
xticks([2:2:24]);
box on;
set(gca,'FontSize',12);
set(gcf,'color','white');


