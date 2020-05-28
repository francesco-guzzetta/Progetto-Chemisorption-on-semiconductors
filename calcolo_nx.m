function [nx, nx_tot] = calcolo_nx(En, psi, kb, T, autovalori, Eg);

for i = 1:autovalori
    f(i) = exp(-(En(i)+(Eg)/2)/(kb*T));
end

[righe, colonne] = size(psi);   
nx = zeros(righe,1);

for i = 1:righe
    for j = 1:colonne
        nx(i) =  nx(i) + f(j)*abs(psi(i,j))^2;                                  
    end                                             
end  

nx_tot = sum(nx);
