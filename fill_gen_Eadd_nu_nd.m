W0=0;
N_sites=12;
epsilon_up(1:N_sites)=0;
epsilon_down=epsilon_up;
x=1:N_sites;
Dt=-0.5:0.05:0.5;
Ec=0:0.25:10;
Nf=1:2*N_sites;
Eg=zeros(length(Ec),length(Dt),length(Nf));
Eadd=zeros(length(Ec),length(Dt),length(Nf));
nu=zeros(length(Ec),length(Dt),length(Nf));
nd=zeros(length(Ec),length(Dt),length(Nf));
for ju=1:length(Ec)
for jt=1:length(Dt)
    fprintf('ju = %d\n',ju);
    fprintf('jt = %d\n',jt);
    dt=Dt(jt);
t_up(1:2:N_sites-1)=-(1-dt);
t_up(2:2:N_sites-1)=-(1+dt);
t_down=t_up;
U(1:N_sites)=Ec(ju);
parfor jn=1:length(Nf)
    [Eg(ju,jt,jn),nu(ju,jt,jn),nd(ju,jt,jn)] = filgs_full(N_sites,Nf(jn),t_up,t_down,epsilon_up,epsilon_down,U,W0,x);
end
Eadd(ju,jt,1)=Eg(ju,jt,1);
for jn=2:length(Nf)
Eadd(ju,jt,jn)=Eg(ju,jt,jn)-Eg(ju,jt,jn-1);
end
end
end
save fill_dat_Eadd_nu_nd