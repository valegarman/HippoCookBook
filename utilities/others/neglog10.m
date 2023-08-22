
function out = neglog10(in,C)
keyboard;
out = sign(in)*(log10(1+abs(in)/(10^C)));

end