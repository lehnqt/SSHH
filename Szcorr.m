function Sc=Szcorr(N_sites,N_up,N_down,w,j,k) %computes the correlation <w|Sz_j Sz_k|w>
w1=(cjck(N_sites,N_up,N_down,w,k,k,1)-cjck(N_sites,N_up,N_down,wg1,k,k,-1))/2;
w2=(cjck(N_sites,N_up,N_down,w1,j,j,1)-cjck(N_sites,N_up,N_down,w1,j,j,-1))/2;
Sc=wg1'*w2;
end