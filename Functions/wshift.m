function w1=wshift(N_sites,N_up,N_down,w,js)
w1=zeros(length(w),1);
basis_up=uperm(de2bi(2^N_sites - 2^(N_sites-N_up)));
if basis_up==0
    basis_up=zeros(1,N_sites);
end

basis_down=uperm(de2bi(2^N_sites - 2^(N_sites-N_down)));
if basis_down==0
    basis_down=zeros(1,N_sites);
end

basis_up_1=basis_up;
basis_down_1=basis_down;

basis_up_1(:,1:js)=circshift(basis_up(:,1:js),[0 1]);
basis_down_1(:,1:js)=circshift(basis_down(:,1:js),[0 1]);

dec_up_1=bi2de(basis_up_1,'left-msb');
dec_down_1=bi2de(basis_down_1,'left-msb');

phase_up=(-1) .^(sum(basis_up(:,1:js-1),2).*basis_up(:,js));
phase_down=(-1) .^(sum(basis_down(:,1:js-1),2).*basis_down(:,js));

size_down=nchoosek(N_sites,N_down);

[~, id_up_1]=sort(dec_up_1);
[~, id_down_1]=sort(dec_down_1);
id_1=(kron(id_up_1,ones(length(id_down_1),1))-1)*size_down+kron(ones(length(id_up_1),1),id_down_1);
phase=kron(phase_up,phase_down);
w=w.*phase;
w1=w(id_1);
end