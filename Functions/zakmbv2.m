function [Eg,flux,zak_phase] = zakmbv2(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,Em)
Nf=400;
Eg=zeros(1,Nf);
flux=linspace(-1,1,Nf)*pi;
H=hbd_dmv2(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,flux(1),Em);
[w1,Eg(1)]=eigs(H,1,'sr');
w_init=w1;
prod=1;
for j=2:length(flux)-1
    H=hbd_dmv2(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,flux(j),Em);
    [w2,Eg(j)]=eigs(H,1,'sr');
   prod=prod*(w1'*w2);
    w1=w2;
end
[sumx_up,sumx_down] = sumx(N_sites,N_up,N_down);
x_op=kron(exp(-1i*(2*pi/N_sites)*sumx_up),exp(-1i*(2*pi/N_sites)*sumx_down));
w2=x_op.*w_init;
prod=prod*(w1'*w2);
zak_phase=angle(prod);
Eg(Nf)=Eg(1);
end