Dt2=-0.4:0.1:0.4;
Sct1=cell(1,length(Dt2));
Sct2=cell(1,length(Dt2));
for k=1:length(Dt2)
    index=find(Dt<Dt2(k),1,'last');
    Sct1{k}=Scd1{length(Ec),index+1};
    Sct2{k}=Scd2{length(Ec),index+1};
end
figs=figure;
for id=1:(length(Dt2)+1)/2
subplot(3,3,id);
colormap('jet');
image(Sct2{id},'CDataMapping','scaled')
set(gca,'Ydir','Normal','FontSize',14);
caxis([-1/4,1/4]);
title(['$\Delta t=$',num2str(Dt2(id))],'interpreter','latex');
set(gcf,'color','white');
end
for id=(length(Dt2)+1)/2+1:length(Dt2)
subplot(3,3,id);
colormap('jet');
image(Sct1{id},'CDataMapping','scaled')
set(gca,'Ydir','Normal','FontSize',14);
caxis([-1/4,1/4]);
title(['$\Delta t=$',num2str(Dt2(id))],'interpreter','latex');
set(gcf,'color','white');
end
set(gcf,'color','white');
figure;
colormap('jet');
caxis([-1/4,1/4]);
colorbar('horiz');
set(gcf,'color','white');
set(gca,'FontSize',14);