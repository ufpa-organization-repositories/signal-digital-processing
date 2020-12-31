% clear all;
% clc;
% close all;

%% parte 1
% carrega os sinais modulados
sinal1_modulado = load('sinal1_Filtrado_normalilzado_dizimado_modulado.mat');
sinal2_modulado = load('sinal2_Filtrado_normalilzado_dizimado_modulado.mat');

% extrai Am dentro do structure
sinal1 = sinal1_modulado.Am;
sinal2 = sinal2_modulado.Am;

% procura o tamanho mínimo
A = length(sinal1_modulado.Am);
B = length(sinal2_modulado.Am);
%[~,tamMin] = min([A B])

% Nao e mais necessario cortar os audios, pois os dois sao do mesmo tamanho
% % corta o audio2 (sinal maior)
% sinal2_modulado.Am(A:B-1) = [];

sinais_somados = sinal1_modulado.Am + sinal2_modulado.Am;
my_fft(sinais_somados, 22050);
title('Espectro dos sinais somados');

% figure();x = fft(sinais_somados); X = abs(x); plot(fftshift(X));

%% Filtro p/ retirar o sinal1 (modulado para 3k Hz)

% fNyquist = 8000/2; % fsampling(que é a freq.max)/2
fNyquist = 22050/2; % fsampling(que é a freq.max)/2
% frequencias de passagem e de corte do filtro passa-faixa
Wp = [2000 4000]/fNyquist; Ws = [1500 4500]/fNyquist; % a faixa(fs1 a fs2) nao pode exceder a faixa filtrada a ser obtida: FALSO
% ripple da banda de passagem e da banda de corte
Rp = 3; Rs = 35;
[n,Wn] = buttord(Wp,Ws,Rp,Rs); % retorna a ordem do vetor e a freq. de corte
[b,a] = butter(n,Wn);
pot_2 = 2^12;%potencia de 2 acima da fNyquist
freqz(b,a,fNyquist,fs2);
sinal1_filtrado_da_soma = filter(b, a, sinais_somados);
my_fft(sinal1_filtrado_da_soma, 8000);title('Resposta em frequencia do sinal filtrado da soma');
% o = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', Ws(1), Wp(1), Wp(2), Ws(2) ,Rs,Rp,Rs);
% od  = design(o,'equiripple');
% freqz(od);
% fvtool(od);

%% Filtro p/ retirar o sinal2 (modulado para 9k Hz)

% fNyquist = 8000/2; % fsampling(que é a freq.max)/2
% frequencias de passagem e de corte do filtro passa-faixa
fNyquist = 22050/2; % fsampling(que é a freq.max)/2
% Wp = [950 1050]/fNyquist; Ws = [900 1100]/fNyquist;
Wp = [8000 10E3]/fNyquist; Ws = [7500 10500]/fNyquist; %OBS: restricao fs2>fNyquist. Na verdade, todas devem ser maiores
% ripple da banda de passagem e da banda de corte
Rp = 3; Rs = 20;
[n,Wn] = buttord(Wp,Ws,Rp,Rs); % retorna a ordem do vetor e a freq. de corte
[b,a] = butter(n,Wn);
pot_2 = 2^12;%potencia de 2 acima da fNyquist
freqz(b,a,fNyquist,fs2);
sinal2_filtrado_da_soma = filter(b, a, sinais_somados);
my_fft(sinal2_filtrado_da_soma, 22050);title('Resposta em frequencia do sinal filtrado da soma');
% o = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', Ws(1), Wp(1), Wp(2), Ws(2) ,Rs,Rp,Rs);
% od  = design(o,'equiripple');
% freqz(od);
% fvtool(od);

%% Parte 2: Demodulacao do sinal 1
% pega parametros para a demodulacao
Ca1_str = load('Ca1.mat');
fa1_str = load('Fa1.mat');
% magn1_str = load('mag1.mat');

% retira os valores da estrutura
Ca1 = Ca1_str.Ca;
fa1 = fa1_str.fa;
% magn1 = magn1_str.magn;

% demodulacao
demod1 = sinal1_filtrado_da_soma.*Ca1;
% demod1 = demod1 - mean(demod1(:));
demod1 = demod1 - 1;
audiowrite('Sinal01_demodulado.wav',demod1, 22050);
figure();x = fft(demod1); X = abs(x); plot(fftshift(X)); title('sinal1_demodulado')

% demod1 = demod1.*magn1;
% audiowrite('Sinal01_demodulado_dc_0_normMagnAmp.wav',demod1, 8000);

%% parte 3: Demodulacao do sinal 2

% pega parametros para a demodulacao
Ca2_str = load('Ca2.mat');
fa2_str = load('Fa2.mat');
magn2_str = load('mag2.mat');

% retira os valores da estrutura
Ca2 = Ca2_str.Ca;
fa2 = fa2_str.fa;
% magn2 = magn2_str.magn;

% procura o tamanho mínimo da portadora
C = length(Ca1);
D = length(Ca2);

%[~,tamMin] = min([A B])

% Nao precisa mais cortar o audio maior pelo tamanho do menor, pois os dois
% tem o mesmo tempo, o mesmo tamanho
% % corta a portadora2 (portadora maior)
% Ca2(C:D-1) = [];

% demodulacao
demod2 = sinal2_filtrado_da_soma.*Ca2;
demod2 = demod2 - 1;
%demod2 = demod2.*magn2;

audiowrite('Sinal02_demodulado.wav',demod2, 22050);
figure();x = fft(demod1); X = abs(x); plot(fftshift(X));
my_fft(demod1, 22050); title('Analise espectral da demodulacao1 ainda em 22050 Hz'); my_fft(demod2, 22050); title('Analise espectral da demodulacao2 ainda em 22050 Hz');
%% Filtro passa-baixa para remover as frequencias adjacente advindas da portadora
% h_str = load('h.mat');
% h = h_str.h;

fp = 2000; %frequ?ncia de passagem
fs = 2500; %frequ?ncia de corte
% fp = 2675; %frequ?ncia de passagem
% fs = 4675; %frequ?ncia de corte

%necess?rio converter para radianos
wp = fp*pi/(22050/2);
ws = fs*pi/(22050/2);

wt = ws - wp; %frequ?ncia de transi??o
% delta_f = abs(ws -wp)/2*pi; %|fs - fp|/fsampling % ta errado
% Comprimento = round(3.3/delta_f); %3.3/wt (Hz) % ta errado
Comprimento = ceil((6.6*pi/wt)) + 1; % tabelado para a janela de hamming
wc = (ws+wp)/2; %frequ?ncia de corte intermedi?ria
hd = my_low_pass_ideal(wc, Comprimento); %fun??o sinc passa baixas ideal
w_hamm = hamming(Comprimento)'; %calcula a janela de hamming
h = hd.*w_hamm; % faz a multiplica??o entre os vetores

freqz(h);
% fvtool(h);


s1df = conv(h, demod1);
s2df = conv(h, demod2);
figure();x = fft(s1df); X = abs(x); plot(fftshift(X)); title('espectro do sinal 1 retirado da soma e eliminada as portadoras');
figure();x = fft(s2df); X = abs(x); plot(fftshift(X)); title('espectro do sinal 2 retirado da soma e eliminada as portadoras');


%% parte 4: expansao do sinal 1 demodulado
% dizimacao por um fator M (Diminuição da taxa de amostragem)
% fsampling = 44100; %freq amost original
% L = 8000/fsampling;
% M = 6; %fs2=8kHz
M = 6;
L = 1/M;

s1exp = s1df(1:M:length(demod1));
audiowrite('sinal01_demodulado_expandido.wav', sinal1_demodulado_expandido, 44100);

%% parte 5: expansao do sinal 2 demodulado
% dizimacao por um fator M (Diminuição da taxa de amostragem)
% fsampling = 44100; %freq amost original
% L = 8000/fsampling;
% M = 6; %fs2=8kHz
M = 6;
L = 1/M;
s2exp = s2df(1:L:length(demod2));
audiowrite('sinal02_demodulado_expandido.wav', sinal2_demodulado_expandido, 44100);