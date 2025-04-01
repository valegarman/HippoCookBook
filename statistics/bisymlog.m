
function y = bisymlog(varargin)
% Bi-symetric blablablalb
% Andrea Gallardo and Manu Valero 2025 
%
x = varargin{1};

if nargin > 1
    C = varargin{2};
else
    C = 0.1;
end

y = sign(x).*(log10(1+abs(x)/(10^C)));

end