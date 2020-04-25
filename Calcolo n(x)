function [nx] = calcolo_nx(En, psi, m0, kT, h, q, Ef, autovalori, ht)

%densita' effettiva di stati nella banda diconduzione
Nc = ((2*pi*m0*kT)/h^2)^1.5; 
A_1D = m0*kT/(pi*ht^2);         




%conversione da elettronVolt a Volt
for i = 1:autovalori
    En(i)=En(i)/q;  
end




%calcolo della quantità di carica nella banda k-esima 
% equazione (10) del paper

for i = 1:autovalori
    Nk(i) = A_1D*Nc*exp(((En(i)-Ef))/kT);                                                                               
end           




% equazione (11) del paper usata per il calcolo della concentrazione totale 
% di elettroni eseguita sommando le concentrazioni delle diverse bande in
% ogni dx così da ricavare una concentrazione variabile nello spazio

[righe, colonne] = size(psi);   
nx = zeros(righe,1);

for i = 1:righe
    for j = 1:colonne
        nx(i) =  nx(i) + Nk(j)*abs(psi(i,j))^2;                                          
    end                                             
end                                                 
                
