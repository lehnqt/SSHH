function [Eg,flux,zak_phase] = spin_zakmbv2(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,Em,sigma)
Nf=100;
Eg=zeros(1,Nf);
flux=linspace(-1,1,Nf)*pi;
mprod=1;
for j=1:length(flux)-1
    H=hbd_dmv2(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,flux(j),Em);
    [w,Eg(j)]=eigs(H,1,'sr');
    rho=spin_partial_trace(N_sites,N_up,N_down,w,sigma);
    mprod=mprod*rho;
end
[sumx_up,sumx_down] = sumx(N_sites,N_up,N_down);
if sigma==1
   x_op=diag(exp(-1i*(2*pi/N_sites)*sumx_up));
elseif sigma==-1
   x_op=diag(exp(-1i*(2*pi/N_sites)*sumx_down));
end
mprod=x_op*mprod;
zak_phase=angle(trace(mprod));
Eg(Nf)=Eg(1);
end