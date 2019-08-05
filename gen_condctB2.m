T=1;
kT=1.38*(10^(-23))*T*10^3/(1.6*10^(-19));
ex_lim=10;
Dmax=6;
chem_min=-20;
chem_max=80;
N_step=200000;  
chem=linspace(chem_min,chem_max,N_step);
EmL=0:0.025:0.5;
g=zeros(length(EmL),length(chem));
paramB2;
Ead=zeros(length(EmL),2*N_sites);
parfor jm=1:length(EmL)
disp(jm);
Em=EmL(jm);
[E0,E_dif,S,N_fill]=filerg(N_sites,t_up,t_down,epsilon_up,epsilon_down,U,W0,x,Em,kappa);
fprintf('calculating transition rate\n');
N_mat=zeros(N_sites+1,N_sites+1);
En=zeros(N_sites+1,N_sites+1,Dmax);
Dn=zeros(N_sites+1,N_sites+1,Dmax);
Q_up=zeros(N_sites+1,N_sites+1,Dmax,Dmax);
Q_down=zeros(N_sites+1,N_sites+1,Dmax,Dmax);
for N_up=0:N_sites
    for N_down=0:N_sites
        N_mat(N_up+1,N_down+1)=N_up+N_down;
        [E,Q1,Q2,D] = transrt_T(N_sites,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,W0,x,Em,kappa,Dmax);
        se=length(E);
        su=size(Q1);
        sd=size(Q2);
        En(N_up+1,N_down+1,1:se)=E;
        Q_up(N_up+1,N_down+1,1:su(1),1:su(2))=Q1;
        Q_down(N_up+1,N_down+1,1:sd(1),1:sd(2))=Q2;
        Dn(N_up+1,N_down+1)=D;
    end
end
fprintf('calculating conductivity\n');
% chem=sort([chem,E_dif]);
gt=zeros(1,length(chem))
for jc=1:length(chem)
    if isempty(find(abs(E_dif-chem(jc))<ex_lim*kT,1))==1
        gt(jc)=0;
        else
    gt(jc)=condct(En,Q_up,Q_down,Dn,N_mat,chem(jc),kT,ex_lim);
    end
end
g(jm,:)=gt;
Ead(jm,:)=E_dif;
end
save dat_condctB_N12_triv
% figcondct=figure;
% plot(chem,g_req);
% title('Conductance ','interpreter','latex')
% xlabel('chemical potential','interpreter','latex') % x-axis label
% ylabel('G','interpreter','latex') % y-axis label
% h=legend(strcat('$\bar{d}=\{',num2str(d1),',',num2str(d2),'\}, \Delta x=',num2str(dx), ', T=',num2str(T),'$'));
% set(h,'interpreter','latex');
% set(gca,'FontSize',18);
% savefig(figcondct,strcat('1dspltc_intersite_Tshape_condct_d_',num2str(d1),'_',num2str(d2),'_dx_',num2str(dx)));
% close(figcondct);
% save(strcat('1dspltc_intersite_Tshape_data_d_',num2str(d1),'_',num2str(d2),'_dx_',num2str(dx)));