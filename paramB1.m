 d2=4.6;
hopp_amp2=6;
d1=6.1;
hopp_amp1=2;
dx=0;
    EB=0;
    V0=0;
    W0=V0;
    EC=40;
    N_sites=12;
    kappa=0;
    
d(1:2:N_sites-1)=d1;
d(2:2:N_sites-1)=d2;
x=zeros(1,N_sites);
for j=2:N_sites
    x(j)=x(j-1)+d(j-1);
end
    x=x+dx*(rand(1,N_sites)-1/2);
for j=1:N_sites-1
    d(j)=x(j+1)-x(j);
end
    epsilon_up(1:N_sites)=-EB;
t_up=zeros(1,N_sites);
t_up(1:2:N_sites-1)=-hopp_amp1;
t_up(2:2:N_sites-1)=-hopp_amp2;
t_down=t_up;
U(1:N_sites)= EC;
for j=1:N_sites
    for k=1:N_sites
        if k~=j
            epsilon_up(j)=epsilon_up(j)-V0/(abs(x(k)-x(j)));
        end
    end
end
epsilon_down=epsilon_up;