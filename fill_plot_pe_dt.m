pen1=zeros(length(Dt),length(Nf));
for jn=1:length(Nf)
    for jt=1:length(Dt)
        pen1(jt,jn)=pen(41,length(Dt)-jt+1,jn);
    end
end
Nf1=ones(length(Dt),1)*Nf;
Dt1=Dt'*ones(1,length(Nf));
Dt2=-0.5:0.001:0.5;
Nf2=ones(length(Dt2),1)*Nf;
Dt2=Dt2'*ones(1,length(Nf));
pen1=interp2(Nf1,Dt1,pen1,Nf2,Dt2);
fig1=figure;
mesh(Nf2,Dt2,pen1);
xlim([0 25]);
colormap(linspecer(100,'sequential'));
box on;
xticks([0,6,12,18,24]);
set(gca,'FontSize',24);
set(gcf,'color','white');
