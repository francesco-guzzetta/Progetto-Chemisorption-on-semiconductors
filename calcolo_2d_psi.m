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

autovalori=10;                   %[/] numero autovalori 

 

 

a = 1e-8;                     %[m] larghezza della buca 

x = linspace(0,a, 1024)';

[qn, fun] = calcolo_2( Eg, kb, T, a, autovalori, m, h,Ef, ht, x);

figure 

plot(qn)

function [qn, fun] = calcolo_2(Eg, kb, T, a,autovalori, m, h,Ef, ht,x)

 

% autovalori esatti

En_esatti = zeros(autovalori, 1);

for i = 1:autovalori

    En_esatti(i)  = (i^2*(h^2))/(8*m*(a^2));

    psi(i,:) = sqrt(2/a)*sin((i*pi*x)/a);

end

 

fun = zeros(autovalori,numel(x));

f = zeros(autovalori,numel(x));

qn = zeros(1,numel(x));

for i = 1:autovalori  

Ei = En_esatti(i);

exp(-(Ei+Ef)/(kb*T))

f(i,:) = kb*T*exp(-(Ei+Ef)./(kb*T)).*(abs(psi(i,:)).^2);

g2d = m/(pi*ht^2);

fun(i,:) = f(i,:).*g2d;

end

figure

plot(x,fun)

 

for i = 1:autovalori

    qn = qn + fun(i,:);

end

end
