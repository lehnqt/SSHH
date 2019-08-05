[E0,E_dif,S,N_fill]=filerg(N_sites,t_up,t_down,epsilon_up,epsilon_down,U,W0,x,Em,kappa);
% figerg=figure;
% plot(N_fill,E_dif,'o');
% title('Addition Energy $E_n-E_{n-1}$','interpreter','latex')
% xlabel('electron number $n$','interpreter','latex') % x-axis label
% ylabel('$\Delta E_n$','interpreter','latex') % y-axis label
% h=legend(strcat('$\bar{d}=\{',num2str(d1),',',num2str(d2),'\}, \Delta x=',num2str(dx),'$'));
% set(h,'interpreter','latex');
% set(gca,'FontSize',18);
% savefig(figerg,strcat('1dspltc_intersite_adderg_d_',num2str(d1),'_',num2str(d2),'_dx_',num2str(dx)));
% close(figerg);