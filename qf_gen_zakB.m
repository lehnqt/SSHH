N_sites=8;
N_up=N_sites/4;
N_down=N_sites/4;
U(1:N_sites)=10;
x=1:N_sites;
flux=0;
Dt=-0.5:0.01:0.5;
epsilon_up(1:N_sites)=0;
epsilon_down=epsilon_up;
Em=0:0.005:0.5;
nug=zeros(length(Em),length(Dt));
ndg=zeros(length(Em),length(Dt));
SzB=zeros(length(Em),length(Dt));
phase=zeros(length(Em),length(Dt));
for jm=1:length(Em)
parfor jt=1:length(Dt)
    fprintf('jm = %d\n',jm);
    fprintf('jt = %d\n',jt);
    dt=Dt(jt);
    [wg,Eg(jm,jt),nug(jm,jt),ndg(jm,jt)]=filgstbcB(N_sites,N_sites/2,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,U,flux,Em(jm));
    SzB(jm,jt)=(nug(jm,jt)-ndg(jm,jt))/2;
    [~,~,phi]=zakmbv2(N_sites,N_up,N_down,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,U,Em(jm));
    phase(jm,jt)=phi/pi;
end
end
save qf_dat_zakB
