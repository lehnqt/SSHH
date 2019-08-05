wnum=input('state number:');
if wnum==1
entt=interp2(entd1,4);
elseif wnum==2
    entt=interp2(entd2,4);
elseif wnum==3
    entt=interp2(entd3,4);
end
colormap(linspecer(100,'sequential'));
% colormap('jet');
image([min(Dt),max(Dt)],[min(Ec),max(Ec)],entt,'CDataMapping','scaled');
set(gca,'Ydir','Normal','FontSize',20);
h=colorbar;
caxis([0 4]);
% xlabel('$\Delta t$','interpreter','latex');
% ylabel('$U$','interpreter','latex');
set(h,'YTick',[0:4]);
box on;
set(gcf,'color','white');