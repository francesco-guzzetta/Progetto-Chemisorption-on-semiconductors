function [En_esatti, psi_esatte] = En_esatti(autovalori, m, a,h, x)
%% autovalori esatti
En_esatti = zeros(autovalori, 1);
for i = 1:autovalori
    En_esatti(i)  = (i^2*(h^2))/(8*m*(a^2));
    psi_esatte(i,:) = sqrt(2/a)*sin((i*pi*x)/a);
end
end