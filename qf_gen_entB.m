load qf_dat_SzB
N_sites=12;
U(1:N_sites)=10;
Dt=-0.5:0.01:0.5;
W0=0;
x=1:N_sites;
epsilon_up(1:N_sites)=0;
epsilon_down=epsilon_up;
N_fill=N_sites/2;
Em=0:0.002:0.2;
entB=zeros(length(Em),length(Dt));
entdB=zeros(length(Em),length(Dt));
for jm=1:length(Em)
parfor jt=1:length(Dt)
    fprintf('jm = %d\n',jm);
    fprintf('jt = %d\n',jt);
    dt=Dt(jt);
N_up=nug(jm,jt);
N_down=ndg(jm,jt);
H=hbd(N_sites,N_up,N_down,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,U,W0,x,Em(jm));
[wg1,Eg]=eigs(H,1,'sa');
    N_sites_L=1;
    entt=0;
   for N_up_L=N_up-min(N_sites-N_sites_L,N_up):min(N_sites_L,N_up)
    for N_down_L=N_down-min(N_sites-N_sites_L,N_down):min(N_sites_L,N_down)
        rho=partial_trace(N_sites,N_up,N_down,N_sites_L,N_up_L,N_down_L,wg1);
        [~,eiv]=eig(rho);
        for j=1:length(eiv)
            if eiv(j,j)~=0
        entt=entt-eiv(j,j)*log2(eiv(j,j));
            end
        end
    end
   end
 entB(jm,jt)=entt;
 
   N_sites_L=2;
   wg1=wshift(N_sites,N_up,N_down,wg1,N_sites);
   entt=0;
for N_up_L=N_up-min(N_sites-N_sites_L,N_up):min(N_sites_L,N_up)
    for N_down_L=N_down-min(N_sites-N_sites_L,N_down):min(N_sites_L,N_down)
        rho=partial_trace(N_sites,N_up,N_down,N_sites_L,N_up_L,N_down_L,wg1);
        [~,eiv]=eig(rho);
        for j=1:length(eiv)
            if eiv(j,j)~=0
        entt=entt-eiv(j,j)*log2(eiv(j,j));
            end
        end
    end
end
entdB(jm,jt)=entt;
end
end
save qf_dat_entB