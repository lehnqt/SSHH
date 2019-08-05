for d_avg=10
for dx=0
    dd=-3;
    d2=d_avg+dd;
    d1=d_avg;
    EB=45.59;
    a=1.7;
    V0=123;
    W0=V0;
    EC=43.86;
    N_sites=10;
    
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
    param;
    adderg;
    reqcondct_plot;
    clearvars -except d_avg dx
end
end