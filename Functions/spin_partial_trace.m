function rho=spin_partial_trace(N_sites,N_up,N_down,w,sigma)
size_up=nchoosek(N_sites,N_up);
size_down=nchoosek(N_sites,N_down);
Mw=reshape(w,[size_down,size_up]);
if sigma==1
    rho=zeros(size_up);
    for j=1:size_up
        for k=j:size_up
            rho(j,k)=sum(Mw(:,j).*conj(Mw(:,k)));
            rho(k,j)=conj(rho(j,k));
        end
    end
    for j=1:size_up
        rho(j,j)=sum(Mw(:,j).*conj(Mw(:,j)));
    end
end    
if sigma==-1
    rho=zeros(size_down);
    for j=1:size_down
        for k=j:size_down
            rho(j,k)=sum(Mw(j,:).*conj(Mw(k,:)));
            rho(k,j)=conj(rho(j,k));
        end
    end
    for j=1:size_down
        rho(j,j)=sum(Mw(j,:).*conj(Mw(j,:)));
    end
end  
end