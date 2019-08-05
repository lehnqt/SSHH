N_sites=12;
U(1:N_sites)=0.01;
x=1:N_sites;
flux=0;
Dt=-0.5:0.01:0.5;
epsilon_up(1:N_sites)=0;
epsilon_down=epsilon_up;
Em=1.5;
nug=zeros(length(Em),length(Dt));
ndg=zeros(length(Em),length(Dt));
SzB=zeros(length(Em),length(Dt));
phase=zeros(length(Em),length(Dt));
for jm=1:length(Em)
parfor jt=1:length(Dt)
    fprintf('jt = %d\n',jt);
    dt=Dt(jt);
    [wg,Eg(jm,jt),nug(jm,jt),ndg(jm,jt)]=filgstbcB(N_sites,N_sites/2,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,U,flux,Em(jm));
    SzB(jm,jt)=(nug(jm,jt)-ndg(jm,jt))/2;
    phase(jm,jt)=zak(N_sites,nug(jm,jt),ndg(jm,jt),wg)/pi;
end
end
plot(Dt,abs(phase));