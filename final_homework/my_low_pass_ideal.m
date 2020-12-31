% Filtro Passa-Baixas ideal
%
%Retorna Sinc(x)
%Para teste de filtro digital passa baixas
function hd = my_low_pass_ideal(wc,M)

alpha = (M-1)/2;
n = (0:M-1);
m = n - alpha + eps;
hd = sin(wc*m) ./(pi*m);

end