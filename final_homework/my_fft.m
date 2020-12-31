function[X, freq] = my_fft(x, Fs)

%x = fft(sinal1_dizimado, 2^17);
%figure();plot(fftshift(abs(x)));

N = length(x);
K = 0:N-1;
T = N/Fs;
freq = K/T;
X = fftn(x)/N;
cutOff = ceil(N/2);
X=X(1:cutOff);
figure();
plot(freq(1:cutOff), abs(X));
title('Espectro de Frequências');
xlabel('Frequencia Hz')
ylabel('Amplitude');

end