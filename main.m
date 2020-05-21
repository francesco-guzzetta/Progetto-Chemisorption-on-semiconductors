clc
clear all
close all
%% MAIN

%%Costanti
q  = 1.602176e-19;  % [C]           carica elettrone
m0 = 9.1095e-31;    % [Kg]          massa elettrone 
h  = 6.6261e-34;    % [J s]         costante di plank 
kb = 1.380645e-23;  % [J K^-1]      costante di boltzman T temperatoura in kelvin
T = 300;            % [k]
ht = h/(2*pi);      % [J s]         costante di plank (tagliata)
mn = 9.11e-31;      % [m]           massa efficace elettrone
Nd = 10e22;         % [m^-3]        concentrazione di drogaggio in un m^3
e0 = 8.85e-12;      % [F m^-1]      epsilon zero
er = 5.4;           % [/]           epsilon r
Eg = 1.08*q;        % [J]
Nc = 2.82e25;      % [m-3]
Nv = 1.04e25;      % [m-3]
Ec = Eg/2 + kb*T*log(sqrt(Nc/Nv)); %[J]
ni = 1.45e16;        %si [m-3]


%% Dati del problema 
a=200e-9;                        %[m] larghezza della buca
autovalori=1;                    %[/] numero autovalori
dx=1.e-11;                       %[m] passo discretizzazione
x = linspace(-a/2,a/2, a/dx)';   %[m] asse x       

%% profilo della buca  di potenziale a pareti infinite 

V = zeros(size(x)); %[V] vettore potenziale
N = length(V);

%% Schrodinger


[En, psi] = Schrodinger_1D(dx, V, autovalori, m0, h, N, x);


%% autovalori esatti
En_esatti = zeros(autovalori, 1);
for i = 1:autovalori
    En_esatti(i)  = (i^2*(h^2))/(8*m0*(a^2));
end


%% calcolo nx

[nx, nx_tot] = calcolo_nx(En, psi, kb, T, autovalori, Ec, Nc);


%% errore

%ritrasformiamo in elettron-volt per maggiore chiarezza 

En = En/q   
En_esatti = En_esatti/q

errore = zeros(autovalori, 1);

for i = 1 : autovalori
    errore(i) = En(i) - En_esatti(i);
end

errore_autov = errore
