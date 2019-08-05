function Dmax=maxdim(N_sites,t_up,t_down,epsilon_up,epsilon_down,U,W0,x,kT,d_lim,ex_lim)
H0=extHbd(N_sites,1,0,t_up,t_down,epsilon_up,epsilon_down,U,W0,x);
D=min(d_lim,length(H0));
E1=eigs(H0,D,'sa');
k=0;
while k<D&&((E1(k+1)-E1(1))/kT <ex_lim)
    k=k+1;
end
D1=k;
H0=extHbd(N_sites,round(N_sites/2),N_sites-round(N_sites/2),t_up,t_down,epsilon_up,epsilon_down,U,W0,x);
D=min(d_lim,length(H0));
E2=eigs(H0,D,'sa');
k=0;
while k<D&&((E2(k+1)-E2(1))/kT <ex_lim)
    k=k+1;
end
D2=k;
if (D1-d_lim)*(D2-d_lim)==0 
    fprintf('dimension limit is too small \n');
end
Dmax=max(D1,D2);
end