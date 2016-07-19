function y=spike(x,L)
%  y=(0.5*(1+cos(pi*min(abs(x)/L,1))));
  y=(1-sin(pi*min(abs(x)/L,1/2)));
end