
function G = condct(En,Q_up,Q_down,Dn,N_mat,chem,kT,ex_lim)
N_sites=length(N_mat)-1;
E_tot=En(:,:,1)-N_mat*chem;
E0=min(E_tot(:));
G=0;
Z=0;
for N_up=0:N_sites
    for N_down=0:N_sites
        if (( En(N_up+1,N_down+1,1)-(N_up+N_down)*chem-E0)/kT)<ex_lim
       for l=1:Dn(N_up+1,N_down+1)
            P=exp(-(1/kT)*(En(N_up+1,N_down+1,l)-(N_up+N_down)*chem-E0));  
            Z=Z+P;
            if N_up>0
           for m=1:Dn(N_up,N_down+1)
               G=G+Q_up(N_up+1,N_down+1,m,l)*P/(1+exp(-(En(N_up+1,N_down+1,l)-En(N_up,N_down+1,m)-chem)/kT));
           end
            end
            if N_down>0
           for m=1:Dn(N_up+1,N_down)
                G=G+Q_down(N_up+1,N_down+1,m,l)*P/(1+exp(-(En(N_up+1,N_down+1,l)-En(N_up+1,N_down,m)-chem)/kT));
           end
            end
       end
        end
    end
end
G=G/(kT*Z);
end
        