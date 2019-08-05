function [wg,Eg,N_up,N_down] = filgstbcB(N_sites,N_fill,t_up,t_down,epsilon_up,epsilon_down,U,flux,Em)
    lower_n=ceil(N_fill/2);
    upper_n=min(N_sites,N_fill);
    for j=lower_n:upper_n
        N_up=j;
        N_down=N_fill-N_up;
      H0=hbd_tbc(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,flux,Em);
      E0=eigs(H0,1,'sa');
        erg(j-lower_n+1)=E0;
    end
[~,id]=min(erg(:));
N_up=id-1+lower_n;
N_down=N_fill-N_up;
[wg, Eg]=eigs(hbd_tbc(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,flux,Em),1,'sa');
end     