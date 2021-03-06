function [En, psi] = Schrodinger_1D(dx, V, autovalori, m0, h, N)

% Matrice H
e = ones(N,1);
c = -(h/2/pi/dx)^2/2/m0;
H = spdiags([c*e -(c*2*e+V) c*e],-1:1, N, N);

%condizione al contorno 
H(1,1)=1;
H(N,N)=1;
H(1,2)=0;
H(N,N-1)=0;

% Calcolo autovalori ed autovettori
[F,D] = eigs(H,autovalori,'smallestabs');   
En = diag(D');                                            

psi = sqrt(1/dx)*F(:,1:autovalori)';  
