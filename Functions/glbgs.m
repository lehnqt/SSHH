function [Eg,N_up,N_down] = glbgs(N_sites,t_up,t_down,epsilon_up,epsilon_down,U,W,chem)
for j=0:N_sites
    for k=0:j
        N_up=j;
        N_down=k;
      H0=extHbd(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,W,chem);
      Eg=eigs(H0,1,'sa');
        erg(j+1,k+1)=Eg;
    end
end
[Eg,L]=min(erg(:));
[N_up, N_down] = ind2sub(size(erg),L);
N_up=N_up-1;
N_down=N_down-1;
end     