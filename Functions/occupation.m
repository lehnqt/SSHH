function [site, density] = occupation(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,W0,x)
H=extHbd(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,W0,x);
[w,~]=eigs(H,1,'sa');
[~,density_up]=partdens(N_sites,N_up,N_down,w,1);
[site,density_down]=partdens(N_sites,N_up,N_down,w,-1);
density=density_up+density_down;
end