function [site, density] = partdens(N_sites,N_up,N_down,w,sigma)
for k=1:N_sites
    site(k)=k;
    density(k)=spdm(N_sites,N_up,N_down,w,k,k,sigma);
end
end

