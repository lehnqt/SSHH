W0=0;
N_sites=12;
epsilon_up(1:N_sites)=0;
epsilon_down=epsilon_up;
x=1:N_sites;
Dt=-0.5:0.05:0.5;
Ec=0:0.25:10;
Nqf=N_sites/2;
lower_n=ceil(Nqf/2);
upper_n=min(N_sites,Nqf);
n1=10;

nE=0;
for k=lower_n:upper_n
        N_up=k;
        N_down=Nqf-N_up;
      szH=nchoosek(N_sites,N_up)*nchoosek(N_sites,N_down);
      nE=nE+min(n1,szH);
end

Eex=zeros(length(Ec),length(Dt),nE);
nuex=zeros(length(Ec),length(Dt),nE);
ndex=zeros(length(Ec),length(Dt),nE);
for ju=1:length(Ec)
    parfor jt=1:length(Dt)
    fprintf('ju = %d\n',ju);
    fprintf('jt = %d\n',jt);
dt=Dt(jt);
id1=1;
id2=0;
nut=zeros(1,nE);
ndt=zeros(1,nE);
Et=zeros(1,nE);
    for k=lower_n:upper_n
        N_up=k;
        N_down=Nqf-N_up;
      szH=nchoosek(N_sites,N_up)*nchoosek(N_sites,N_down);
      H0=extHbd(N_sites,N_up,N_down,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,Ec(ju)*ones(1,N_sites),W0,x);
      neig=min(n1,szH);
      id2=id2+neig;
      E0=eigs(H0,neig,'sa');
      nut(id1:id2)=N_up*ones(1,neig);
      ndt(id1:id2)=N_down*ones(1,neig);
      Et(id1:id2)=E0;
     id1=id2+1;
    end
[m, im]=sort(Et);    
Eex(ju,jt,:)=m;
nuex(ju,jt,:)=nut(im);
ndex(ju,jt,:)=ndt(im);
    end
end
save qf_dat_Eex
