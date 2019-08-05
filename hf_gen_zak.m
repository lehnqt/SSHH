N_sites=12;
x=1:N_sites;
flux=0;
Em=0;
Dt=-0.5:0.01:0.5;
Ec=0:0.1:10;
epsilon_up(1:N_sites)=0;
epsilon_down=epsilon_up;
nug=zeros(length(Ec),length(Dt));
ndg=zeros(length(Ec),length(Dt));
SzB=zeros(length(Ec),length(Dt));
phase=zeros(length(Ec),length(Dt));
for ju=1:length(Ec)
parfor jt=1:length(Dt)
    fprintf('ju = %d\n',ju);
    fprintf('jt = %d\n',jt);
    dt=Dt(jt);
    U=Ec(ju)*ones(1,N_sites);
    [wg,Eg(ju,jt),nug(ju,jt),ndg(ju,jt)]=filgstbcB(N_sites,N_sites,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,U,flux,Em);
    SzB(ju,jt)=(nug(ju,jt)-ndg(ju,jt))/2;
    phase(ju,jt)=zak(N_sites,nug(ju,jt),ndg(ju,jt),wg)/pi;
end
end
save hf_dat_zak
