function [En,Q_up,Q_down,Dn] = transrt_T(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,W0,x,Em,kappa,Dmax)
    Hn=hbd(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,W0,x,Em,kappa);
    Dn=min(length(Hn),Dmax);
    [wn,En]=eigs(Hn,Dn,'sa');
    En=diag(En);
    if N_up==0
        Q_up=zeros(1,Dn);
    else
    Hnp=hbd(N_sites,N_up-1,N_down,t_up,t_down,epsilon_up,epsilon_down,U,W0,x,Em,kappa);
    Dnp=min(length(Hnp),Dmax);
    [wnp,~]=eigs(Hnp,Dnp,'sr');
    for j=1:Dnp
        for k=1:Dn
            vec=0;
            for p=1:1
            vec=vec+cj(N_sites,N_up,N_down,wn(:,k),p,1);
            end
            Q_up(j,k)=((wnp(:,j))'*vec)^2;
        end
    end
    end
    if N_down==0
        Q_down=zeros(1,Dn);
    else
    Hnp=hbd(N_sites,N_up,N_down-1,t_up,t_down,epsilon_up,epsilon_down,U,W0,x,Em,kappa);
    Dnp=min(length(Hnp),Dmax);
    [wnp,~]=eigs(Hnp,Dnp,'sa');
    for j=1:Dnp
        for k=1:Dn
            vec=0;
            for p=1:1
            vec=vec+cj(N_sites,N_up,N_down,wn(:,k),p,-1);
            end
            Q_down(j,k)=((wnp(:,j))'*vec)^2;
        end
    end
    end
end
    