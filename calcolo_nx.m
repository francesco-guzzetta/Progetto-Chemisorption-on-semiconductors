function [nx_tot] = calcolo_nx (En, Ef, kb, T);

  f = exp (-(En-Ef)/(kb*T));
  nx_tot = sum (f);

end
