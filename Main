clc
clear all
%% MAIN

%Costanti
q  = 1.68e-19;      % [C]           carica elettrone
m0 = 9.1095e-31;    % [Kg]          massa elettrone 
h  = 6.6261e-34;    % [J s]         costante di plank 
kT = 25.8e-3;       % [J K^-1 K]    costante di boltzman T temperatoura in kelvin
ht = h/(2*pi);      % [J s]         costante di plank (tagliata)
mn = 1.91e-31;      % [m]           massa efficace elettrone
Ef = -0.58;         % [V]           potenziale di fermi
e0 = 8.85e-12;      % [F m^-1]      epsilon zero
er = 5.4;           % [/]           epsilon r
Nd = 10e22;         % [cm^-3]       concentrazione di drogaggio in un cm^3


%% Dati del problema 
a=5e-9;             %[m] larghezza della buca
autovalori=5;       %[1] numero autovalori
dx=5.e-11;          %[m] passo discretizzazione
x = -a/2:dx:a/2';   %[m] asse x

%% profilo della buca  di potenziale a pareti infinite 


V = zeros(size(x)); %[V] vettore potenziale
N = length(V);

%% Schrodinger

[En, psi] = Schrodinger_1D(dx, V, autovalori, m0, h, N);


%% autovalori esatti
En_esatti = zeros(autovalori, 1);
for i = 1:autovalori
    En_esatti(i)  = (i^2*(h^2))/(8*m0*(a^2));
end

%% calcolo nx

% calcolo nx
[nx] = calcolo_nx(En, psi, m0, kT, h, q, Ef, autovalori, ht);


%calcolo analitico (soluzione esatta) 
[nx_esatto] = calcolo_nx_esatto(En_esatti, psi, m0, kT, h, Ef, autovalori, ht, q);




%% errore

%ritrasformiamo in elettron-volt per maggiore chiarezza 
En = En/q   
En_esatti = En_esatti/q

errore = zeros(autovalori, 1);

for i = 1 : autovalori
    errore(i) = En(i) - En_esatti(i);
end

errore
