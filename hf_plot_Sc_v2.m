Dt2=[-0.48,0,0.5];
Sct=cell(1,length(Dt2));
for k=1:length(Dt2)
    index=find(Dt<Dt2(k),1,'last');
    Sct{k}=Scr1{length(Ec),index+1};
end
Sct{3}=(Scr1{length(Ec),length(Dt)}+Scr2{length(Ec),length(Dt)}+2*Scr3{length(Ec),length(Dt)})/4;
figs=figure;
for id=1:length(Dt2)
figure;
colormap('jet');
image(Sct{id},'CDataMapping','scaled')
set(gca,'Ydir','Normal','FontSize',18);
h=colorbar();
caxis([-1/4,1/4]);
set(h,'ytick',[-1/4,0,1/4]);
% title(['$\Delta t=$',num2str(Dt2(id))],'interpreter','latex');
set(gcf,'color','white');
end