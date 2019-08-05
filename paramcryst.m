epsilon_up(1:N_sites)=-EB;
t_up=zeros(1:N_sites);
t_up(1:2:N_sites-1)=-hopp_amp1;
t_up(2:2:N_sites)=-hopp_amp2;
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