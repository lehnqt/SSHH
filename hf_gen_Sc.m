W0=0;
N_sites=12;
epsilon_up(1:N_sites)=0;
epsilon_down=epsilon_up;
x=1:N_sites;
Dt=-0.5:0.05:0.5;
Ec=0:0.25:10;
N_up=N_sites/2;
N_down=N_sites/2;
Scr1=cell(length(Ec),length(Dt));
Scr2=cell(length(Ec),length(Dt));
Scr3=cell(length(Ec),length(Dt));
for ju=1:length(Ec)
    parfor jt=1:length(Dt)
            fprintf('ju = %d\n',ju);
            fprintf('jt = %d\n',jt);
dt=Dt(jt);
H=extHbd(N_sites,N_up,N_down,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,Ec(ju)*ones(1,N_sites),W0,x);
[wg,Eg]=eigs(H,2,'sa');
wg1=wg(:,1);
wg2=wg(:,2);
H=extHbd(N_sites,N_up+1,N_down-1,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,Ec(ju)*ones(1,N_sites),W0,x);
[wg,Eg3]=eigs(H,1,'sa');
wg3=wg;
Scrt1=zeros(N_sites,N_sites);
Scrt2=zeros(N_sites,N_sites);
Scrt3=zeros(N_sites,N_sites);
for j=1:N_sites
    for k=1:N_sites
    w1=(cjck(N_sites,N_up,N_down,wg1,k,k,1)-cjck(N_sites,N_up,N_down,wg1,k,k,-1))/2;
    w2=(cjck(N_sites,N_up,N_down,w1,j,j,1)-cjck(N_sites,N_up,N_down,w1,j,j,-1))/2;
    Scrt1(j,k)=wg1'*w2;
    w1=(cjck(N_sites,N_up,N_down,wg2,k,k,1)-cjck(N_sites,N_up,N_down,wg2,k,k,-1))/2;
    w2=(cjck(N_sites,N_up,N_down,w1,j,j,1)-cjck(N_sites,N_up,N_down,w1,j,j,-1))/2;
    Scrt2(j,k)=wg2'*w2;
    w1=(cjck(N_sites,N_up+1,N_down-1,wg3,k,k,1)-cjck(N_sites,N_up+1,N_down-1,wg3,k,k,-1))/2;
    w2=(cjck(N_sites,N_up+1,N_down-1,w1,j,j,1)-cjck(N_sites,N_up+1,N_down-1,w1,j,j,-1))/2;
    Scrt3(j,k)=wg3'*w2;
    end
end
Scr1{ju,jt}=Scrt1;
Scr2{ju,jt}=Scrt2;
Scr3{ju,jt}=Scrt3;
    end
end
save hf_dat_Sc
            