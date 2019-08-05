N_sites=12;
U(1:N_sites)=10;
Dt=-0.5:0.01:0.5;
W0=0;
x=1:N_sites;
epsilon_up(1:N_sites)=0;
epsilon_down=epsilon_up;
N_fill=N_sites/2;
Em=0:0.002:0.2;
Eg=zeros(length(Em),length(Dt));
nug=zeros(length(Em),length(Dt));
ndg=zeros(length(Em),length(Dt));
SzB=zeros(length(Em),length(Dt));
for jm=1:length(Em)
parfor jt=1:length(Dt)
    fprintf('jm = %d\n',jm);
    fprintf('jt = %d\n',jt);
    dt=Dt(jt);
[Eg(jm,jt),nug(jm,jt),ndg(jm,jt)]=filgsB(N_sites,N_sites/2,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,U,W0,x,Em(jm));
SzB(jm,jt)=(nug(jm,jt)-ndg(jm,jt))/2;
end
end
save qf_dat_SzB