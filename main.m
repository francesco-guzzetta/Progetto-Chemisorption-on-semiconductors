clc
clear all
close all 

%% MAIN
 
%Costanti
q  = 1.60217662e-19;    % [C]           carica elettrone
mh = 1.6398e-31;        % [Kg]          massa lacuna 
h  = 6.6261e-34;        % [J s]         costante di plank 
kb = 1.38064852e-23;    % [J K^-1]      costante di boltzman 
T = 300;                % [k]           T temperatoura in kelvin
ht = 1.0545718e-34;     % [J s]         costante di plank (tagliata)
m0 = 9.9359e-31;        % [m]           massa elettrone
Nd = 1e22;              % [m^-3]        concentrazione di drogaggio in un m^3
e0 = 8.85e-12;          % [F m^-1]      costante dielettrica nel vuoto
er = 5.4;               % [/]           
Eg = 1.1082*q;          % [J]
Nc = sqrt (((m0*kb*T/(ht^2*pi))^3) / 2);        % [m-3]
Nv = sqrt (((mh*kb*T/(ht^2*pi))^3) / 2);        % [m-3]
ni = sqrt(Nc * Nv) * exp (-(Eg/2)/(kb*T));      % [m-3]

%% Se NON drogato

Ef0 = Eg/2+kb*T*log(sqrt(Nc/Nv)); %[J]   %si [m-3]

%% Se Drogato
%Ef = Ef0 - kb*T*log((Nd/(ni)))

% dati buca 
buca = logspace(-9, -7, 20);

ni = sqrt(Nc * Nv) * exp (-(Eg/2)/(kb*T))
for i = 1 : numel(buca) 
a = buca(i);                     %[m] larghezza della buca                               
dx=3.e-12;                       %[m] passo discretizzazione
x = linspace(0,a, a/dx)';   %[m] asse x   
V = zeros(size(x));              %[V] vettore potenziale
N = length(V);
autovalori = 100;
Ef = Ef0;

%En esatti
%[En, psi] = En_esatti(autovalori, m, a, h, x);

%En schrodinger
[En, psi] = Schrodinger_1D(dx, V, autovalori, m0, h, N);

[qn, fun] = calcolo_n(Eg, kb, T, a,autovalori, m0,Ef, ht, x, En ,psi);

hold on
xlabel('x/L')
ylabel('Densità di portatori')
plot(x/a,qn);
plot([0 1],[ni ni])
end 
