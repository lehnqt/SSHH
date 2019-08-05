density=zeros(2*N_sites,N_sites);
for j=1:2*N_sites
    N_up=ceil(j/2);
    N_down=j-N_up;
    [site,density(j,:)]=occupation(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,W0,x);
end