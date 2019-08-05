N_sites_L=1;
entrp=0;
for N_up_L=N_up-min(N_sites-N_sites_L,N_up):min(N_sites_L,N_up)
    for N_down_L=N_down-min(N_sites-N_sites_L,N_down):min(N_sites_L,N_down)
        rho=partial_trace(N_sites,N_up,N_down,N_sites_L,N_up_L,N_down_L,w);
        [~,eiv]=eig(rho);
        for j=1:length(eiv)
            if eiv(j,j)~=0
        entrp=entrp-eiv(j,j)*log2(eiv(j,j));
            end
        end
    end
end