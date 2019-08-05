function psi_j = cdagj(N_sites, N_up, N_down, w, j,sigma)
size_up=nchoosek(N_sites,N_up);
size_down=nchoosek(N_sites,N_down);
if sigma==1
    basis_up=uperm(de2bi(2^N_sites - 2^(N_sites-N_up),N_sites));
    basis_up_2=uperm(de2bi(2^N_sites - 2^(N_sites-N_up-1),N_sites));
    dec_up=bi2de(basis_up,'left-msb');
    dec_up_2=bi2de(basis_up_2,'left-msb');
    size_up_2=nchoosek(N_sites,N_up+1);
    psi_j=zeros(size_up_2*size_down,1);
   
    [~,row_index_up,column_index_up]=intersect(dec_up + 2^(N_sites-j)*(1-basis_up(:,j)),dec_up_2);
    exponent = (-1) .^(sum(basis_up_2(column_index_up,1:j),2));
    
   for j_down=1:size_down
   row_index=(row_index_up-1)*size_down + j_down;
   column_index=(column_index_up-1)*size_down+j_down;
   psi_j(column_index)= w(row_index) .* exponent;
   end
   
elseif sigma== -1
    basis_down=uperm(de2bi(2^N_sites - 2^(N_sites-N_down),N_sites));
    basis_down_2=uperm(de2bi(2^N_sites - 2^(N_sites-N_down-1),N_sites));
    dec_down=bi2de(basis_down,'left-msb');
    dec_down_2=bi2de(basis_down_2,'left-msb');
    size_down_2=nchoosek(N_sites,N_down+1);
    psi_j=zeros(size_down_2*size_up,1);
   
    [~,row_index_down,column_index_down]=intersect(dec_down + 2^(N_sites-j)*(1-basis_down(:,j)),dec_down_2);
    exponent = (-1) .^(sum(basis_down_2(column_index_down,1:j),2));
    
    
    for j_up =1:size_up
   row_index=(j_up-1)*size_down + row_index_down;
   column_index=(j_up-1)*size_down_2+column_index_down;
    psi_j(column_index)= w(row_index) .* exponent;
   end
   
else
    error('The last variable must be -1 for spin down or 1 for spin up.');
end

