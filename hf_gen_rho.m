W0=0;
N_sites=12;
epsilon_up(1:N_sites)=0;
epsilon_down=epsilon_up;
x=1:N_sites;
N_up=ceil(N_sites/2);
N_down=floor(N_sites/2);
dt=0.5;
Ec=10;
H=extHbd(N_sites,N_up,N_down,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,Ec*ones(1,N_sites),W0,x);
[wg,Eg]=eigs(H,2,'sa');
wg1=wg(:,1);
wg2=wg(:,2);
H=extHbd(N_sites,N_up+1,N_down-1,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,Ec*ones(1,N_sites),W0,x);
[wg3,~]=eigs(H,1,'sa');
N_sites_L=2;
wed1=wshift(N_sites,N_up,N_down,wg1,N_sites);
wed2=wshift(N_sites,N_up,N_down,wg2,N_sites);
wed1=wshift(N_sites,N_up,N_down,wed1,2);
wed2=wshift(N_sites,N_up,N_down,wed2,2);
rhoed1=cell(3,3);
rhoed2=cell(3,3);
rhoed3=cell(3,3);
for N_up_L=N_up-min(N_sites-N_sites_L,N_up):min(N_sites_L,N_up)
    for N_down_L=N_down-min(N_sites-N_sites_L,N_down):min(N_sites_L,N_down)
        rhoed1{N_up_L+1,N_down_L+1}=partial_trace(N_sites,N_up,N_down,N_sites_L,N_up_L,N_down_L,wed1);
        rhoed2{N_up_L+1,N_down_L+1}=partial_trace(N_sites,N_up,N_down,N_sites_L,N_up_L,N_down_L,wed2);
    end
end
wed3=wshift(N_sites,N_up+1,N_down-1,wg3,N_sites);
for N_up_L=N_up+1-min(N_sites-N_sites_L,N_up+1):min(N_sites_L,N_up+1)
    for N_down_L=N_down-1-min(N_sites-N_sites_L,N_down-1):min(N_sites_L,N_down-1)
        rhoed3{N_up_L+1,N_down_L+1}=partial_trace(N_sites,N_up+1,N_down-1,N_sites_L,N_up_L,N_down_L,wed3);
    end
end
rhodm1=cell(N_sites/2-1,3,3);
rhodm2=cell(N_sites/2-1,3,3);
rhodm3=cell(N_sites/2-1,3,3);
for k=1:(N_sites/2-1)
    js=2*k+1;
wdm1=wshift(N_sites,N_up,N_down,wg1,js);
wdm1=wshift(N_sites,N_up,N_down,wdm1,js);
wdm2=wshift(N_sites,N_up,N_down,wg2,js);
wdm2=wshift(N_sites,N_up,N_down,wdm2,js);
for N_up_L=N_up-min(N_sites-N_sites_L,N_up):min(N_sites_L,N_up)
    for N_down_L=N_down-min(N_sites-N_sites_L,N_down):min(N_sites_L,N_down)
        rhodm1{k,N_up_L+1,N_down_L+1}=partial_trace(N_sites,N_up,N_down,N_sites_L,N_up_L,N_down_L,wdm1);
        rhodm2{k,N_up_L+1,N_down_L+1}=partial_trace(N_sites,N_up,N_down,N_sites_L,N_up_L,N_down_L,wdm2);
    end
end

wdm3=wshift(N_sites,N_up+1,N_down-1,wg3,js);
wdm3=wshift(N_sites,N_up+1,N_down-1,wdm3,js);
for N_up_L=N_up+1-min(N_sites-N_sites_L,N_up+1):min(N_sites_L,N_up+1)
    for N_down_L=N_down-1-min(N_sites-N_sites_L,N_down-1):min(N_sites_L,N_down-1)
        rhodm3{k,N_up_L+1,N_down_L+1}=partial_trace(N_sites,N_up+1,N_down-1,N_sites_L,N_up_L,N_down_L,wdm3);
    end
end
end

dt=-0.5;
H=extHbd(N_sites,N_up,N_down,-(1+dt*(-1).^x),-(1+dt*(-1).^x),epsilon_up,epsilon_down,Ec*ones(1,N_sites),W0,x);
[wgt,Eg]=eigs(H,1,'sa');
rhodm1t=cell(N_sites/2,3,3);
for k=1:(N_sites/2)
    js=2*k;
wdm1t=wshift(N_sites,N_up,N_down,wgt,js);
wdm1t=wshift(N_sites,N_up,N_down,wdm1t,js);
for N_up_L=N_up-min(N_sites-N_sites_L,N_up):min(N_sites_L,N_up)
    for N_down_L=N_down-min(N_sites-N_sites_L,N_down):min(N_sites_L,N_down)
        rhodm1t{k,N_up_L+1,N_down_L+1}=partial_trace(N_sites,N_up,N_down,N_sites_L,N_up_L,N_down_L,wdm1t);
    end
end
end

save dat_hf_rho