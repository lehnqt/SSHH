function psi_jk = cjck(N_sites, N_up, N_down, w, j, k,sigma) %cdagjck
size_up=nchoosek(N_sites,N_up);
size_down=nchoosek(N_sites,N_down);
l1=min(j,k);
l2=max(j,k);
if sigma==1
    basis_up=uperm(de2bi(2^N_sites - 2^(N_sites-N_up)));
    dec_up=bi2de(basis_up,'left-msb');
    psi_jk=zeros(length(w),1);
   
    if N_up>0
    [~,row_index_up,column_index_up]=intersect((dec_up+2^(N_sites-j)- 2^(N_sites-k)).*(1-basis_up(:,j)*(1-eq(j,k))).* basis_up(:,k),dec_up);
    exponent = (-1) .^(sum(basis_up(column_index_up,l1:l2),2)-1);
    
   for j_down=1:size_down
   row_index=(row_index_up-1)*size_down + j_down;
   column_index=(column_index_up-1)*size_down+j_down;
   psi_jk(column_index)= w(row_index) .* exponent;
   end
    end
elseif sigma== -1
    basis_down=uperm(de2bi(2^N_sites - 2^(N_sites-N_down)));
    dec_down=bi2de(basis_down,'left-msb');
        psi_jk=zeros(length(w),1);
        
        if N_down>0
    [~,row_index_down,column_index_down]=intersect((dec_down+2^(N_sites-j)- 2^(N_sites-k)).*(1-basis_down(:,j)*(1-eq(j,k))).* basis_down(:,k),dec_down);
    exponent = (-1) .^(sum(basis_down(column_index_down,l1:l2),2)-1);

   for j_up =1:size_up
   row_index=(j_up-1)*size_down + row_index_down;
   column_index=(j_up-1)*size_down+column_index_down;
    psi_jk(column_index)= w(row_index) .* exponent;
   end
        end
else
    error('The last variable must be -1 for spin down or 1 for spin up.');
end

