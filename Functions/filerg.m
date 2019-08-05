function [E0,E_dif,S,N_fill]=filerg(N_sites,t_up,t_down,epsilon_up,epsilon_down,U,W0,x,Em,kappa)
for k=1:2*N_sites
N_fill(k)=k;
[E0(k),N_up,N_down]=filgs(N_sites,N_fill(k),t_up,t_down,epsilon_up,epsilon_down,U,W0,x,Em,kappa);
S(k)=(N_up-N_down)/2;
end
E_dif(1)=E0(1);
for k=2:2*N_sites
E_dif(k)=E0(k)-E0(k-1);
end
end
 %scatter(N_fill,E0)
%title('Ground state energy','interpreter','latex')
%xlabel('electron number','interpreter','latex') % x-axis label
%ylabel('Energy','interpreter','latex') % y-axis label
%ax = gca;
%ax.XTick = (0:1:12);
%gtext('$N_{sites}=6,t= 1, U=10$','interpreter','latex')
