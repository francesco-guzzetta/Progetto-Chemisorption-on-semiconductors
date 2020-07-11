function [qn, fun] = calcolo_n(Eg, kb, T, a,autovalori, m, Ef, ht, x, En ,psi)

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