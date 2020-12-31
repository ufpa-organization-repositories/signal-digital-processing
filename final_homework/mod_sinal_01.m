clear all;
clc;
close all;
% [audio1, fsampling] = audioread('Sinal_E.wav');
% [audio1, fsampling] = audioread('nota_C5.wav');
[audio1, fsampling] = audioread('boate_azul.wav');


%%%
% fp = 3000; %frequ?ncia de passagem
% fs = 4350; %frequ?ncia de corte
fp = 4000; %frequ?ncia de passagem
fs = 18050; %frequ?ncia de corte
% fp = 2675; %frequ?ncia de passagem
% fs = 4675; %frequ?ncia de corte

%necess?rio converter para radianos
wp = fp*pi/(fsampling/2);
ws = fs*pi/(fsampling/2);

wt = ws - wp; %frequ?ncia de transi??o
% delta_f = abs(ws -wp)/2*pi; %|fs - fp|/fsampling % ta errado
% Comprimento = round(3.3/delta_f); %3.3/wt (Hz) % ta errado
Comprimento = ceil((6.6*pi/wt)) + 1; % tabelado para a janela de hamming
wc = (ws+wp)/2; %frequ?ncia de corte intermedi?ria
hd = my_low_pass_ideal(wc, Comprimento); %fun??o sinc passa baixas ideal
w_hamm = hamming(Comprimento)'; %calcula a janela de hamming
h = hd.*w_hamm; % faz a multiplica??o entre os vetores

freqz(h);
fvtool(h);

sinal_filtrado = conv(h, audio1); %convolu??o entre o filtro e o sinal
% sinal_filtrado = sinal_filtrado*1E3;
audiowrite('Sinal01_filtrado.wav',sinal_filtrado, fsampling);
a =fft(audio1); A = abs(a); plot(fftshift(A)); figure();fa = fft(sinal_filtrado); FA = abs(fa); plot(fftshift(FA))
my_fft(audio1, 44100); my_fft(sinal_filtrado, 44100);

% %normalização do sinal
% g = abs(sinal_filtrado); %aplicando modulo no sinal filtrado
% magn = max(g); %passando para a variavel magn a magnitude maxima do sinal.
% norm = sinal_filtrado./magn; %realizando o processo de normalização
% z = norm+1; %sinal "norm" normalizado + 1, isso faz com que o sinal tenha um deslocamento na amplitude,
% ma=max(abs(norm));
% sinal_filtrado=z;
% audiowrite('Sinal01_filtrado_normalizado.wav',sinal_filtrado, fsampling);

sinal_filtrado = sinal_filtrado + 1;

% dizimacao por um fator M (Diminuição da taxa de amostragem)
% freq_nyquist = fsampling/2;
% fs2 = freq_nyquist; %nova freq amost
fs2 = 22050; %nova freq amostragem
M = round(fsampling/fs2);

sinal1_dizimado = sinal_filtrado(1:M:length(sinal_filtrado));
audiowrite('Sinal01_dizimado.wav',sinal1_dizimado, fs2); % tira do nivel DC 1 e volta pro nivel DC 0. Agora podemos dar o ganho

% modulacao
% t = 0.0001:0.0001:8.6611;
% t = 0.0001:0.0001:25.9832;
% t = 0:1/fsampling:1.96395;
% t = 0:1/8000:9.1876; % equivale a linha de codigo abaixo
aux = [0:length(sinal1_dizimado)-1];
t = aux./fs2;
% t = linspace(0,length(sinal1_dizimado)/fsampling,length(sinal1_dizimado));
fa = 3E3;
A = 1;

Ca = A*sin(2*pi*fa*t);
Am = (sinal1_dizimado').*Ca;

figure(1)
subplot(3,1,1);
plot(t,Ca)
title('Portadora do Sinal 01')
subplot(3,1,2)
plot(t,sinal1_dizimado)
title('Sinal modulante (Sinal 01)')
subplot(3,1,3)
plot(t,Am)
title('Sinal modulado (Sinal 01)');

my_fft(Ca, fs2);
title('Frequencia da portadora do Sinal 01');
%figure(); plot(angle(Ca)); title('Fase da portadora');

my_fft(sinal1_dizimado, fs2);
title('Frequencia do sinal modulante (Sinal 01)');

my_fft(Am, fs2);
title('Frequencia do sinal modulado (Sinal 01)');
audiowrite('sinal1_Filtrado_normalilzado_dizimado_modulado.wav',Am, fs2);

save('sinal1_Filtrado_normalilzado_dizimado_modulado.mat','Am');
save('Ca1.mat','Ca');
save('Fa1.mat','fa');
% save('mag1.mat','magn');
save('h.mat', 'h');
% 
% %Demodulacao
% demod = Am .* Ca;
% 
% my_fft(demod, fsampling)
% title('Frequencia do sinal demodulado (Sinal 01)');
% audiowrite('Sinal01_demodulado.wav',demod,fsampling);
% 
% figure(); plot(Xs); title('sinal modulante (Sinal 01)'); figure(); plot(demod); title('sinal demodulado (Sinal 01)')

