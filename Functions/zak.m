function zak_phase = zak(N_sites,N_up,N_down,wg)
N_dimer=N_sites/2;
[sumx_up,sumx_down] = sumx_dimer(N_sites,N_up,N_down);
x_op=kron(exp(1i*(2*pi/N_dimer)*sumx_up),exp(1i*(2*pi/N_dimer)*sumx_down));
w=x_op.*wg;
prod=wg'*w;
zak_phase=angle(prod);
end