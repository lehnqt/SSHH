function [sumx_up,sumx_down] = sumx(N_sites,N_up,N_down)

basis_up=uperm(de2bi(2^N_sites - 2^(N_sites-N_up)));
if basis_up==0
    basis_up=zeros(1,N_sites);
end
sumx_up=zeros(length(basis_up(:,1)),1);
for j=1:N_sites
    sumx_up=sumx_up+basis_up(:,j)*j;
end

%similarly for spin down
basis_down=uperm(de2bi(2^N_sites - 2^(N_sites-N_down)));
if basis_down==0
    basis_down=zeros(1,N_sites);
end
sumx_down=zeros(length(basis_down(:,1)),1);
for j=1:N_sites
    sumx_down=sumx_down+basis_down(:,j)*j;
end

end