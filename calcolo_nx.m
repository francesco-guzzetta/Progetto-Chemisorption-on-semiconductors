function [nx, nx_tot] = calcolo_nx(En, psi, kb, T, autovalori, Ec, Nc)

%% 1D

%abbiamo utilizzato il valore di Ec come valore di riferimento così che
%al crescere della larghezza della banda il primo autovalore effettivamente
%converge a zero e quindi il valore di nx viene coerente

for i = 1:autovalori
    Nk(i) = Nc*exp((-(Ec+En(i)))/(kb*T));
end  

[righe, colonne] = size(psi);   
nx = zeros(righe,1);

for i = 1:righe
    for j = 1:colonne
        nx(i) =  nx(i) + Nk(j)*abs(psi(i,j))^2;                                  
    end                                             
end  

nx_tot = sum(nx)
