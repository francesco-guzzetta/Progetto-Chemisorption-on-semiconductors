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


%% Se NON drogato

Ef0 = Eg/2+kb*T*log(sqrt(Nc/Nv)); %[J]   %si [m-3]

% dati buca 
buca = logspace(-9, -7, 20);

%buca = 1e-7
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
[En, psi] = Schrodinger_1D(dx, V, autovalori, m, h, N, x, a);

[qn, fun] = calcolo_2(Eg, kb, T, a,autovalori, m,Ef, ht, x, En ,psi);

hold on
xlabel('x/L')
plot(x/a,qn);
plot([0 1],[ni ni])

end 




function [qn, fun] = calcolo_2(Eg, kb, T, a,autovalori, m, Ef, ht, x, En ,psi)

fun = zeros(autovalori,numel(x));
f = zeros(autovalori,numel(x));
qn = zeros(1,numel(x));

for i = 1:autovalori  
Ei = En(i);
f(i,:) = kb*T*exp(-(Ei+Ef)./(kb*T)).*(abs(psi(i,:)).^2);
g2d = m/(pi*ht^2);
fun(i,:) = f(i,:).*g2d;
end

for i = 1:autovalori
    qn = qn + fun(i,:);
end
end


function [En_esatti, psi_esatte] = En_esatti(autovalori, m, a,h, x)
%% autovalori esatti
En_esatti = zeros(autovalori, 1);
for i = 1:autovalori
    En_esatti(i)  = (i^2*(h^2))/(8*m*(a^2));
    psi_esatte(i,:) = sqrt(2/a)*sin((i*pi*x)/a);
end
end
