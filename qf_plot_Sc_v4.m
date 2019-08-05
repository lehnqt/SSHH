Dt2=[-0.48,0,0.5];
Sct1=cell(1,length(Dt2));
Sct2=cell(1,length(Dt2));
for k=1:length(Dt2)
    index=find(Dt<Dt2(k),1,'last');
    Sct1{k}=Scd1{length(Ec),index+1};
    Sct2{k}=Scd2{length(Ec),index+1};
end
figs=figure;
for id=1:(length(Dt2)+1)/2
figure;
colormap('jet');
image(Sct2{id},'CDataMapping','scaled')
set(gca,'Ydir','Normal','FontSize',18);
caxis([-1/4,1/4]);
title(['$\Delta t=$',num2str(Dt2(id))],'interpreter','latex');
set(gcf,'color','white');
h=colorbar();
set(h,'ytick',[-1/4,0,1/4]);
end
for id=(length(Dt2)+1)/2+1:length(Dt2)
figure;
colormap('jet');
image(Sct1{id},'CDataMapping','scaled')
set(gca,'Ydir','Normal','FontSize',18);
caxis([-1/4,1/4]);
title(['$\Delta t=$',num2str(Dt2(id))],'interpreter','latex');
set(gcf,'color','white');
h=colorbar();
set(h,'ytick',[-1/4,0,1/4]);
end