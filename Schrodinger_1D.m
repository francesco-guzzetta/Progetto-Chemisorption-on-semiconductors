function [En, psi] = Schrodinger_1D(dx, V, autovalori, m0, h, N, x)

% Costruzione della matrice per la risoluzione dell'equazione quadratica
% con il metodo delle differenze finite con codice ottimizzato

%% Matrice H
e = ones(N,1);
H = -(h/2/pi/dx)^2/2/m0*spdiags([1*e -2*e 1*e],-1:1, N, N);

%condizione al contorno 
H(1,1)=1;
H(N,N)=1;

% Calcolo autovalori ed autovettori
[F,D] = eigs(H,autovalori,'smallestabs');   %restituisce la matrice diagonale D con gli autovalori e
                                            %restituisce una matrice F le cui colonne corrisondono
                                            %agli autovettori (destri) tali che H*V = V*D
             
W = diag(D);     % W è un vettore contenente gli autovalori    
                
[En, kk] = sort(W,'ascend');  %viene restituito il vettore E che è 
                              %il vettore W in ordine crescente
                              %kk è il vettore delle posizioni degli
                              %elementi,non ordinate,nel vettore w
                              %quindi E=W(kk)
                              
En = En(1:autovalori);
psi_completa = F(:,kk);       %ordina le colonne della matrice F in base
                              %al vettore di ordinamento kk

psi = psi_completa(:,1:autovalori); %si prendono della matrice psi_completa solo 
                                    %le "autofunzioni" richieste(dal valore autovalori)
                                    


