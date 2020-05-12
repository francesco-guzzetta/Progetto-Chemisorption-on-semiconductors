 function [nx, En_tot] = calcolo_nx(En, psi, m0, kb, T, h, Ef, autovalori, ht, q)

%densita' effettiva di stati nella banda diconduzione
Nc = ((2*pi*m0*kb*T)/h^2)^1.5;  %[(1/m-3)^3/2]     % 1.2547e+25
A_1D = m0*kb*T/(pi*ht^2);       %                  % 1.0799e+17



%% 3d

[righe, colonne] = size(psi);
En_tot = 0;
count = 1;
psi_tot = zeros(righe, 1);
En_tot_temp = 0;

for i = 1:autovalori
    for j = 1:autovalori
        for k = 1:autovalori
            En_tot_temp = En(i)+En(j)+En(k);
            En_tot_temp = round(En_tot_temp, 27);
            if max(En_tot == En_tot_temp) == 1
                count = count;
            else 
                En_tot(count) = En_tot_temp;
                psi_tot(:,count) = psi(:,i).*psi(:,j).*psi(:,k);
                count = count +1;
            end
        end
    end
end

[En_tot, kk] = sort(En_tot,'ascend')
psi_tot = psi_tot(:,kk);


% calcolo della quantità di carica nella banda k-esima 
% equazione (10) del paper

for i = 1:length(En_tot)
    Nk(i) = 2.82*10^25*exp(-((En_tot(i)-Ef))/(kb*T));  
end  

% equazione (11) del paper usata per il calcolo della concentrazione totale 
% di elettroni eseguita sommando le concentrazioni delle diverse bande in
% ogni dx così da ricavare una concentrazione variabile nello spazio

[righe, colonne] = size(psi_tot);   
nx = zeros(righe,1);

for i = 1:righe
    for j = 1:colonne
        nx(i) =  nx(i) + Nk(j)*abs(psi_tot(i,j))^2;                                          
    end                                             
end            




%% 1D
% 
% for i = 1:autovalori
%     Nk(i) = 2.82*10^25*exp(-((En(i)-Ef))/(kb*T));  
%     %Nk(i) = A_1D*Nc*exp(-((En(i)-Ef))/(kb*T));
% end  
% 
% % equazione (11) del paper usata per il calcolo della concentrazione totale 
% % di elettroni eseguita sommando le concentrazioni delle diverse bande in
% % ogni dx così da ricavare una concentrazione variabile nello spazio
% 
% [righe, colonne] = size(psi);   
% nx = zeros(righe,1);
% 
% for i = 1:righe
%     for j = 1:colonne
%         nx(i) =  nx(i) + Nk(j)*abs(psi(i,j))^2;                                          
%     end                                             
% end                                                 
% 

 
% for i = 1:autovalori
%     En(i)=En(i);  
% end
