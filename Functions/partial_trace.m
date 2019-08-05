function rho=partial_trace(N_sites,N_up,N_down,N_sites_L,N_up_L,N_down_L,w)
N_sites_R=N_sites-N_sites_L;
N_up_R=N_up-N_up_L;
N_down_R=N_down-N_down_L;

if N_up_L<=min(N_sites_L,N_up)&&N_up_L>=N_up-min(N_sites_R,N_up)&&N_down_L<=min(N_sites_L,N_down)&&N_down_L>=N_down-min(N_sites_R,N_down)
basis_up=uperm(de2bi(2^N_sites - 2^(N_sites-N_up)));
if basis_up==0
    basis_up=zeros(1,N_sites);
end
dec_up=bi2de(basis_up,'left-msb');
basis_down=uperm(de2bi(2^N_sites - 2^(N_sites-N_down)));
if basis_down==0
    basis_down=zeros(1,N_sites);
end
dec_down=bi2de(basis_down,'left-msb');
 basis_up_R=uperm(de2bi(2^(N_sites_R) - 2^(N_sites_R-N_up_R)));
if basis_up_R==0
    basis_up_R=zeros(1,N_sites_R);
end
dec_up_R=bi2de(basis_up_R,'left-msb');

basis_down_R=uperm(de2bi(2^(N_sites_R) - 2^(N_sites_R-N_down_R)));
if basis_down_R==0
    basis_down_R=zeros(1,N_sites_R);
end
dec_down_R=bi2de(basis_down_R,'left-msb');

dec_up_sub=bi2de(basis_up(:,N_sites_L+1:N_sites),'left-msb');
dec_down_sub=bi2de(basis_down(:,N_sites_L+1:N_sites),'left-msb');

size_down=nchoosek(N_sites,N_down);
size_rho=nchoosek(N_sites_L,N_up_L)*nchoosek(N_sites_L,N_down_L);
rho=zeros(size_rho);
for ku=1:size(dec_up_R)
    for kd=1:size(dec_down_R)
id_up=find(dec_up_sub==dec_up_R(ku));
id_down=find(dec_down_sub==dec_down_R(kd));
idw=(kron(id_up,ones(length(id_down),1))-1)*size_down+kron(ones(length(id_up),1),id_down);
wA=w(idw);
rho=rho+wA*wA';
    end
end
else
    rho=0;
end
end