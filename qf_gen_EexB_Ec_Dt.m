N_sites=12;
Ec=0:0.5:10;
Dt=-0.5:0.05:0.5;
W0=0;
x=1:N_sites;
epsilon_up(1:N_sites)=0;
epsilon_down=epsilon_up;
N_fill=N_sites/2;
lower_n=ceil(N_fill/2);
upper_n=min(N_sites,N_fill);
nmax=10;
nE=0;
Em=0:0.01:1.1;
fprintf('Em max = %d\n',Em(length(Em)));
for k=lower_n:upper_n
        N_up=k;
        N_down=N_fill-N_up;
      szH=nchoosek(N_sites,N_up)*nchoosek(N_sites,N_down);
      nE=nE+min(nmax,szH);
end
EexB=zeros(length(Em),length(Dt),length(Ec),nE);
nuexB=zeros(length(Em),length(Dt),length(Ec),nE);
ndexB=zeros(length(Em),length(Dt),length(Ec),nE);
for jm=1:length(Em)
    fprintf('jm = %d\n',jm);
    for jt=1:length(Dt)
         fprintf('jt = %d\n',jt);
parfor jc=1:length(Ec)
    fprintf('jc = %d\n',jc);
    U=Ec(jc)*ones(1,N_sites);
    dt=Dt(jt);
    Eextt=zeros(nE,1);
    nut=zeros(nE,1);
    ndt=zeros(nE,1);
    id1=1;
    id2=0;
    for k=lower_n:upper_n
        N_up=k;
        N_down=N_fill-N_up;
      szH=nchoosek(N_sites,N_up)*nchoosek(N_sites,N_down);
      H0=hbd(N_sites,N_up,N_down,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,U,W0,x,Em(jm));
      neig=min(nmax,szH);
      id2=id2+neig;
      E0=eigs(H0,neig,'sa');
      nut(id1:id2)=N_up*ones(1,neig);
      ndt(id1:id2)=N_down*ones(1,neig);
      Eextt(id1:id2)=E0;
     id1=id2+1;
    end
[m, im]=sort(Eextt);    
EexB(jm,jt,jc,:)=m;
nuexB(jm,jt,jc,:)=nut(im);
ndexB(jm,jt,jc,:)=ndt(im);
end
    end
end
save qf_dat_EexB_Ec_Dt