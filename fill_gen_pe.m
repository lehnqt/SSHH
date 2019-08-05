load fill_dat_Eadd_nu_nd
Nmu=100;
mu=zeros(length(Ec),length(Dt),Nmu);
pen=zeros(length(Ec),length(Dt),length(Nf));
pemu=zeros(length(Ec),length(Dt),Nmu);
ptotmu=zeros(length(Ec),length(Dt),Nmu);
for ju=1:length(Ec)
for jt=1:length(Dt)
    fprintf('ju = %d\n',ju);
    fprintf('jt = %d\n',jt);
    dt=Dt(jt);
t_up(1:2:N_sites-1)=-(1-dt);
t_up(2:2:N_sites-1)=-(1+dt);
t_down=t_up;
U(1:N_sites)=Ec(ju);
Et=zeros(1,length(Nf));
parfor jn=1:length(Nf)
    Et(jn)=Eadd(ju,jt,jn);
    H=extHbd(N_sites,nu(ju,jt,jn),nd(ju,jt,jn),t_up,t_down,epsilon_up,epsilon_down,U,W0,x); 
    [wf,erg]=eigs(H,1,'sa');
   pen(ju,jt,jn)=spdm(N_sites,nu(ju,jt,jn),nd(ju,jt,jn),wf,1,1,1)+spdm(N_sites,nu(ju,jt,jn),nd(ju,jt,jn),wf,1,1,-1)+spdm(N_sites,nu(ju,jt,jn),nd(ju,jt,jn),wf,N_sites,N_sites,1)+spdm(N_sites,nu(ju,jt,jn),nd(ju,jt,jn),wf,N_sites,N_sites,-1);
end
    mut=linspace(min(Et)-3,max(Et)+3,Nmu-length(Et));
    mut=sort([mut,Et]);
    mu(ju,jt,:)=mut;
for jm=1:length(mut)
    jn=find(Et<mut(jm),1,'last');
     if ~isempty(jn)
pemu(ju,jt,jm)=pen(ju,jt,jn);
ptotmu(ju,jt,jm)=Nf(jn);
     end
end
end
end
save fill_dat_pe
