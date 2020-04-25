function [En, psi] = Schrodinger_1D(dx, V, autovalori, m0, h, N)

% Costruzione della matrice per la risoluzione dell'equazione quadratica
% con il metodo delle differenze finite

H = diag([(-(h/2/pi/dx)^2/2/m0)*ones(1,N-2) 0],-1) + ...
    diag([0 (-(h/2/pi/dx)^2/2/m0)*ones(1,N-2)],1) + ...
    diag(((h/2/pi/dx)^2/m0+V));

%condizione al contorno 
H(1,1)=1;
H(N,N)=1;

% Calcolo autovalori ed autovettori
[F,D] = eig(H);  %restituisce la matrice diagonale D con gli autovalori e
                 %restituisce una matrice F le cui colonne corrisondono
                 %agli autovettori (destri) tali che H*V = V*D
             
W = diag(D);     % W è un vettore contenente gli autovalori    
                
[En, kk] = sort(W,'ascend');  %viene restituito il vettore E che è 
                              %il vettore W in ordine crescente
                              %kk è il vettore delle posizioni degli
                              %elementi,non ordinate,nel vettore w
                              %quindi E=W(kk)
                              
En = En(1:autovalori);        %seleziono solo gli autovalori richiesti

psi_completa = F(:,kk);       %ordina le colonne della matrice F in base
                              %al vettore di ordinamento kk

psi = psi_completa(:,1:autovalori); %si prendono della matrice psi_completa solo 
                                    %le "autofunzioni" richieste(dal valore autovalori)
                                    
%% grafici

% un primo esempio delle prime 5 autofunzioni per avere un'idea 
%(da migliorare con label e correzzione della grandezza del dominio) 

figure                                    
subplot(2, 3, 1), plot(psi(:,1))
subplot(2, 3, 2), plot(psi(:,2))
subplot(2, 3, 3), plot(psi(:,3))
subplot(2, 3, 4), plot(psi(:,4))
subplot(2, 3, 5), plot(psi(:,5))



