clc
clear all
close all 
%% MAIN

%Costanti
q  = 1.60217662e-19;    % [C]           carica elettrone
mh = 1.6398e-31;        % [Kg]          massa elettrone 
h  = 6.6261e-34;        % [J s]         costante di plank 
kb = 1.38064852e-23;    % [J K^-1]      costante di boltzman 
T = 300;                % [k]           T temperatoura in kelvin
ht = 1.0545718e-34;     % [J s]         costante di plank (tagliata)
m = 9.9359e-31;         % [m]           massa efficace elettrone
Nd = 1e22;              % [m^-3]        concentrazione di drogaggio in un m^3
e0 = 8.85e-12;          % [F m^-1]      epsilon zero
er = 5.4;               % [/]           epsilon r
Eg = 1.1082*q;          % [J]
Nc = sqrt (((m*kb*T/(ht^2*pi))^3) / 2);         % [m-3]
Nv = sqrt (((mh*kb*T/(ht^2*pi))^3) / 2);        % [m-3]
ni = 1.45*10^16;                                % [m-3]


%% Se NON drogato
Ef = Eg/2+kb*T*log(sqrt(Nc/Nv)); %[J]   %si [m-3]
 
%% Dati del problema 
autovalori=100;                   %[/] numero autovalori 

%% Dati del problema 
autovalori=100;                   %[/] numero autovalori 

buca = logspace(-9, -7, 20);
for j = 1 : numel(buca) 
a = buca(j);                     %[m] larghezza della buca 

[qn(j)] = calcolo_2(Eg, kb, T, a,autovalori, m, h,Ef, ht);
end
qn
semilogy(qn,'LineWidth',3)

function [qn] = calcolo_2(Eg, kb, T, a,autovalori, m, h,Ef, ht)

% autovalori esatti
En_esatti = zeros(autovalori, 1);
for i = 1:autovalori
    En_esatti(i)  = (i^2*(h^2))/(8*m*(a^2));
end
fun = 0;
q1=zeros(1,autovalori);
for i = 1:autovalori  
Ei = En_esatti(i);
f =  kb*T*exp(-(Ei+Ef)./(kb*T));
g2d = m/(pi*ht^2);
fun = fun + g2d.*f;
end
qn=fun/a^2;
end
