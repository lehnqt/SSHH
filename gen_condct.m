warning('off','all')
T=1;
kT=1.38*(10^(-23))*T*10^3/(1.6*10^(-19));
d_lim=5;
ex_lim=5;
Dmax=5;
gen_adderg;

chem_min=min(E_dif)-10;
chem_max=max(E_dif)+10;
fprintf('calculating transition rate\n');
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
N_step=100000;  
chem=linspace(chem_min,chem_max,N_step);
% chem=sort([chem,E_dif]);
g=zeros(1,length(chem));
for k=1:length(chem)
    if isempty(find(abs(E_dif-chem(k))<ex_lim*kT,1))==1
        g(k)=0;
        else
    g(k)=condct(En,Q_up,Q_down,Dn,N_mat,chem(k),kT,ex_lim);
    end
end

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