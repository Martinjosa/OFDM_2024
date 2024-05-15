% Programa general del sistema UWB-OFDM-SAR%
clc;
clear all;
close all;
%======================== Inicio de variables ============================%
numBits= 100000;                                  % 100 mil bits
SNR=[1:30];
for k=1:20
    for i=1:length(SNR)
        %i
        x= randi([0,1],numBits,1);                        % Genero el vector columna de 100 mil b
        bpsk_signal = pskmod(x,2);                        % Modulacion bpsk
        bpsk_awgn=awgn(bpsk_signal,SNR(i),'measured');
        bpsk_r=pskdemod(bpsk_awgn,2);
        isErr=x~=bpsk_r;
        BER(i,k)= nnz(isErr)/numBits;
    end
end
Ber_promedio=mean(BER,2);
figure;
title("Rendimiento en canal awgn")
xlabel("SNR[dB]")
ylabel("Bit Error Rate")
semilogy (SNR, Ber_promedio)
grid on
