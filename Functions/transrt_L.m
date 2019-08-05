function [En,Q_up,Q_down,Dn] = transrt_L(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,W0,x,Dmax)
    Hn=extHbd(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,W0,x);
    Dn=min(length(Hn),Dmax);
    [wn,En]=eigs(Hn,Dn,'sa');
    En=diag(En);
    if N_up==0
        Q_up=zeros(1,Dn);
    else
    Hnp=extHbd(N_sites,N_up-1,N_down,t_up,t_down,epsilon_up,epsilon_down,U,W0,x);
    Dnp=min(length(Hnp),Dmax);
    [wnp,~]=eigs(Hnp,Dnp,'sa');
    for j=1:Dnp
        for k=1:Dn
            Ql=((wnp(:,j))'*cj(N_sites,N_up,N_down,wn(:,k),1,1))^2;
            Qr=((wnp(:,j))'*cj(N_sites,N_up,N_down,wn(:,k),N_sites,1))^2;
            Q_up(j,k)=2/(1/Ql+1/Qr);
        end
    end
    end
    if N_down==0
        Q_down=zeros(1,Dn);
    else
    Hnp=extHbd(N_sites,N_up,N_down-1,t_up,t_down,epsilon_up,epsilon_down,U,W0,x);
    Dnp=min(length(Hnp),Dmax);
    [wnp,~]=eigs(Hnp,Dnp,'sa');
    for j=1:Dnp
        for k=1:Dn
            Ql=((wnp(:,j))'*cj(N_sites,N_up,N_down,wn(:,k),1,-1))^2;
            Qr=((wnp(:,j))'*cj(N_sites,N_up,N_down,wn(:,k),N_sites,-1))^2;
            Q_down(j,k)=2/(1/Ql+1/Qr);
        end
    end
    end
end
    