clear all;
clc;
close all;
% [audio2, fsampling] = audioread('Sinal_F.wav');
% [audio2, fsampling] = audioread('nota_A5.wav');
[audio2, fsampling] = audioread('tocando_em_frente.wav');

%%%
% fp = 3000; %frequ?ncia de passagem
% fs = 4350; %frequ?ncia de corte
% fp = 2675; %frequ?ncia de passagem
% fs = 4675; %frequ?ncia de corte
fp = 4000; %frequ?ncia de passagem
fs = 18050; %frequ?ncia de corte

%necess?rio converter para radianos
wp = fp*pi/(fsampling/2);
ws = fs*pi/(fsampling/2);

wt = ws - wp; %frequ?ncia de transi??o
% delta_f = abs(ws -wp)/2*pi; %|fs - fp|/fsampling
% Comprimento = round(3.3/delta_f);
Comprimento = ceil((6.6*pi/wt)) + 1; %3.3/wt (Hz)
wc = (ws+wp)/2; %frequ?ncia de corte intermedi?ria
hd = my_low_pass_ideal(wc, Comprimento); %fun??o sinc passa baixas ideal
w_hamm = hamming(Comprimento)'; %calcula a janela de hamming
h = hd.*w_hamm; % faz a multiplica??o entre os vetores

freqz(h);
fvtool(h);

sinal_filtrado = conv(h, audio2); %convolu??o entre o filtro e o sinal
% sinal_filtrado = sinal_filtrado*1E3;
audiowrite('Sinal02_filtrado.wav',sinal_filtrado, fsampling);

% %normalização do sinal
% g = abs(sinal_filtrado); %aplicando modulo no sinal filtrado
% magn = max(g); %passando para a variavel magn a magnitude maxima do sinal.
% norm = sinal_filtrado./magn; %realizando o processo de normalização
% z = norm + 1; %sinal "norm" normalizado + 1, isso faz com que o sinal tenha um deslocamento na amplitude,
% %ma=max(abs(norm));
% sinal_filtrado=z;
% audiowrite('Sinal02_filtrado_normalizado.wav',sinal_filtrado, fsampling);

sinal_filtrado = sinal_filtrado + 1;

% dizimacao por um fator M (Diminuição da taxa de amostragem)
% freq_nyquist = fsampling/2;
% fs2 = freq_nyquist; %nova freq amost
% M = round(fsampling/fs2);
fs2 = 22050; %nova freq amostragem
M = round(fsampling/fs2);


sinal2_dizimado = sinal_filtrado(1:M:length(sinal_filtrado));
audiowrite('Sinal02_dizimado.wav',sinal2_dizimado, fs2); % sinal2_dizimado - 1 para voltar para o nivel dc 0

%modulacao

% t = 0.0001:0.0001:9.9155;
% t = 0.0001:0.0001:29.7464;
% t = 0:1/8000:9.1876; % ajustado de acordo com o sinal 2 que tem duracao de pouco menos de 13 segundos
aux = [0:length(sinal2_dizimado)-1];
t = aux./fs2;
% t = linspace(0,length(sinal2_dizimado)/fsampling,length(sinal2_dizimado));
fa = 9E3;
A = 1;

Ca = A*sin(2*pi*fa*t);
Am = (sinal2_dizimado').*Ca;

figure(1)
subplot(3,1,1);
plot(t,Ca)
title('Portadora do Sinal 02')
subplot(3,1,2)
plot(t,sinal2_dizimado)
title('Sinal modulante (Sinal 02)')
subplot(3,1,3)
plot(t,Am)
title('Sinal modulado (Sinal 02)');

my_fft(Ca, fs2);
title('Frequencia da portadora do Sinal 02');
%figure(); plot(angle(Ca)); title('Fase da portadora');

my_fft(sinal2_dizimado, fs2);
title('Frequencia do sinal modulante (Sinal 02)');

my_fft(Am, fs2);
title('Frequencia do sinal modulado (Sinal 02)');

audiowrite('sinal2_Filtrado_normalilzado_dizimado_modulado.wav',Am, fs2);

save('sinal2_Filtrado_normalilzado_dizimado_modulado.mat','Am');

save('Ca2.mat','Ca');
save('Fa2.mat','fa');
% save('mag2.mat','magn');

% %Demodulacao
% demod = Am .* Ca;
% 
% my_fft(demod, fsampling)
% title('Frequencia do sinal demodulado (Sinal 01)');
% audiowrite('Sinal01_demodulado.wav',demod,fsampling);
% 
% figure(); plot(Xs); title('sinal modulante (Sinal 01)'); figure(); plot(demod); title('sinal demodulado (Sinal 01)')

