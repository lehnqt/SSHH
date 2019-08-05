function H_tot = hbd(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,W0,x,Em,kappa)
%Version 4 fix the bug in hopp_amplitude and has a more transparent code.
%double count when N_sites=2, valid only for N_sites>2
kappa=kappa/N_sites;
epsilon_up=epsilon_up-Em*ones(1,N_sites);
epsilon_down=epsilon_down+Em*ones(1,N_sites);
t_up(1:N_sites-1)=t_up(1:N_sites-1)*exp(-1i*kappa);
t_up(N_sites)=t_up(N_sites)*exp(1i*kappa);
t_down(1:N_sites-1)=t_down(1:N_sites-1)*exp(-1i*kappa);
t_down(N_sites)=t_down(N_sites)*exp(1i*kappa);

basis_up=uperm(de2bi(2^N_sites - 2^(N_sites-N_up)));
if basis_up==0
    basis_up=zeros(1,N_sites);
end

dec_up=bi2de(basis_up,'left-msb');
row_index_up=[];
column_index_up=[];
hopp_amplitude_up=[];
if N_up>0
j=1;
%c+_j c_j+1 shift, t is the hopping in front of  c_jc+_j+1
 while j < N_sites
    [~,L1,L2]=intersect((dec_up+2^(N_sites-j)-2^(N_sites-j-1)).*(1-basis_up(:,j)).* basis_up(:,j+1),dec_up);
row_index_up=[row_index_up;L1];
column_index_up=[column_index_up;L2];
hopp_amplitude_up= [hopp_amplitude_up; ones(length(L1),1)*t_up(j)];
j=j+1;
 end
 % t_Nsites is the hopping in front of  c_1c+_N
[~,L1,L2]=intersect((dec_up+2^(N_sites-1)-2^0).*(1-basis_up(:,1)).* basis_up(:,N_sites),dec_up);
row_index_up=[row_index_up;L1];
column_index_up=[column_index_up;L2];
exponent_up = (-1) .^(sum(basis_up(L2,1:N_sites),2)-1);
hopp_amplitude_up= [hopp_amplitude_up; exponent_up .* ones(length(L1),1)*t_up(N_sites)];
end

%similarly for spin down
basis_down=uperm(de2bi(2^N_sites - 2^(N_sites-N_down)));
if basis_down==0
    basis_down=zeros(1,N_sites);
end

dec_down=bi2de(basis_down,'left-msb');
row_index_down=[];
column_index_down=[];
hopp_amplitude_down=[];
if N_down>0
j=1;
while j < N_sites
     [~,L1,L2]=intersect((dec_down+2^(N_sites-j)-2^(N_sites-j-1)).*(1-basis_down(:,j)).* basis_down(:,j+1),dec_down);
row_index_down=[row_index_down;L1];
column_index_down=[column_index_down;L2];
hopp_amplitude_down= [hopp_amplitude_down; ones(length(L1),1)*t_down(j)];
j=j+1;
end
[~,L1,L2]=intersect((dec_down+2^(N_sites-1)-2^0).*(1-basis_down(:,1)).* basis_down(:,N_sites),dec_down);
row_index_down=[row_index_down;L1];
column_index_down=[column_index_down;L2];
exponent_down = (-1) .^(sum(basis_down(L2,1:N_sites),2)-1);
hopp_amplitude_down= [hopp_amplitude_down; exponent_down .* ones(length(L1),1)*t_down(N_sites)];
end

%Generate the kinetic part of the Hamiltonian
size_up=nchoosek(N_sites,N_up);
size_down=nchoosek(N_sites,N_down);
H_up=sparse([row_index_up; column_index_up], [column_index_up; row_index_up],[hopp_amplitude_up; conj(hopp_amplitude_up)],size_up,size_up);
H_down=sparse([row_index_down; column_index_down], [column_index_down; row_index_down],[hopp_amplitude_down; conj(hopp_amplitude_down)],size_down,size_down);
I_up=speye(size_up);
I_down=speye(size_down);
T1=kron(H_up,I_down);
T2=kron(I_up,H_down);
H_tot=T1+T2;
nj_up=cell(1,N_sites);
nj_down=cell(1,N_sites);
for j=1:N_sites
    nj_up{j}=spdiags(basis_up(:,j),0,size_up,size_up);
    nj_up_tensor=kron(epsilon_up(j)*nj_up{j},I_down);
    nj_down{j}=spdiags(basis_down(:,j),0,size_down,size_down);
    nj_down_tensor=kron(I_up,epsilon_down(j)*nj_down{j});
    Uj=kron(U(j)*nj_up{j},nj_down{j});
    H_tot=H_tot+nj_up_tensor+nj_down_tensor+Uj;
end
if W0~=0
 for k=1:N_sites-1
     for j=k+1:N_sites
         W=W0/abs(x(j)-x(k));
                W1=kron(W*nj_up{j},nj_down{k});
                 W2=kron(W*nj_up{k},nj_down{j});
                  W3=kron(nj_up{j}*nj_up{k},W*I_down);
                   W4=kron(W*I_up,nj_down{j}*nj_down{k});
    H_tot=H_tot+sparse(W1)+sparse(W2)+sparse(W3)+sparse(W4);
     end
 end
end
end