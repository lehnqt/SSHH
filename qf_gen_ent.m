W0=0;
N_sites=12;
epsilon_up(1:N_sites)=0;
epsilon_down=epsilon_up;
x=1:N_sites;
Dt=-0.5:0.05:0.5;
Ec=0:0.25:10;
ent=zeros(length(Ec),length(Dt));
entd=zeros(length(Ec),length(Dt));
N_up=ceil(N_sites/4);
N_down=floor(N_sites/4);
for ju=1:length(Ec)
 parfor jt=1:length(Dt)
    fprintf('ju = %d\n',ju);
    fprintf('jt = %d\n',jt);
    dt=Dt(jt);
    H=extHbd(N_sites,N_up,N_down,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,Ec(ju)*ones(1,N_sites),W0,x);
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
 ent(ju,jt)=entt;
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
entd(ju,jt)=entt;
end
end
save qf_dat_ent

 