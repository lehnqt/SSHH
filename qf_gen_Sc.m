W0=0;
N_sites=12;
epsilon_up(1:N_sites)=0;
epsilon_down=epsilon_up;
x=1:N_sites;
Dt=-0.5:0.05:0.5;
Ec=0:0.25:10;
N_up=N_sites/4;
N_down=N_sites/4;
Scr=cell(length(Ec),length(Dt));
Scd1=cell(length(Ec),length(Dt));
Scd2=cell(length(Ec),length(Dt));
for ju=1:length(Ec)
    parfor jt=1:length(Dt)
            fprintf('ju = %d\n',ju);
            fprintf('jt = %d\n',jt);
dt=Dt(jt);
H=extHbd(N_sites,N_up,N_down,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,Ec(ju)*ones(1,N_sites),W0,x);
[wg,Eg]=eigs(H,1,'sa');
Scrt=zeros(N_sites,N_sites);
Scdt1=zeros(N_sites/2+1,N_sites/2+1);
Scdt2=zeros(N_sites/2,N_sites/2);
for j=1:N_sites
    for k=1:N_sites
    w1=(cjck(N_sites,N_up,N_down,wg1,k,k,1)-cjck(N_sites,N_up,N_down,wg1,k,k,-1))/2;
    w2=(cjck(N_sites,N_up,N_down,w1,j,j,1)-cjck(N_sites,N_up,N_down,w1,j,j,-1))/2;
    Scrt(j,k)=wg1'*w2;
    end
end
Scdt1(1,1)=Scrt(1,1);
Scdt1(N_sites/2+1,N_sites/2+1)=Scrt(N_sites,N_sites);
Scdt1(1,N_sites/2+1)=Scrt(1,N_sites);
Scdt1(N_sites/2+1,1)=Scrt(N_sites,1);
for k=2:(N_sites/2)
        Scdt1(1,k)=Scrt(1,2*k-2)+Scrt(1,2*k-1);
        Scdt1(k,1)=Scrt(2*k-2,1)+Scrt(2*k-1,1);
        Scdt1(k,N_sites/2+1)=Scrt(2*k-2,N_sites)+Scrt(2*k-1,N_sites);
        Scdt1(N_sites/2+1,k)=Scrt(N_sites,2*k-2)+Scrt(N_sites,2*k-1);
end
for j=2:N_sites/2
    for k=2:N_sites/2
   Scdt1(j,k)=Scrt(2*j-2,2*k-2)+Scrt(2*j-2,2*k-1)+Scrt(2*j-1,2*k-2)+Scrt(2*j-1,2*k-1);
    end
end
for j=1:N_sites/2
    for k=1:N_sites/2
   Scdt2(j,k)=Scrt(2*j-1,2*k-1)+Scrt(2*j-1,2*k)+Scrt(2*j,2*k-1)+Scrt(2*j,2*k);
    end
end
Scr{ju,jt}=Scrt;
Scd1{ju,jt}=Scdt1;
Scd2{ju,jt}=Scdt2;
    end
end
save qf_dat_Sc
            