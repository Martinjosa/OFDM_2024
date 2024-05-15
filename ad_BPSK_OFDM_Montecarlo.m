% Programa general del sistema UWB-OFDM-SAR%
clc;
clear all;
close all;
%======================== Inicio de variables ============================%
numBits= 100000;                                  % 100 mil bits
SNR=[1:30];
num_subportadoras = 64; % Número de subportadoras en OFDM
cyclic_prefix_length = 16; % Longitud del prefijo cíclico en OFDM
Mod = 2; % Bits por símbolo [1 2 4]

% Ajustar el número de bits para que sea un múltiplo del número de subportadoras
num_bits_s = ceil(numBits / num_subportadoras) * num_subportadoras;

for k=1:20
    for i=1:length(SNR)
        x= randi([0,1],num_bits_s,1);                % Genero el vector columna, bits en paralelo
        bpsk_signal = pskmod(x,Mod);                   % Modulacion bpsk
        
        % Generar símbolos OFDM
        num_symbols = ceil(num_bits_s / num_subportadoras); % Número de simbolos OFDM necesarios
        pulsos_ofdm = reshape(bpsk_signal, num_subportadoras, num_symbols); % Agrupar los bits en símbolos OFDM (las 64 portadoras en paralelo)
        
        % Modulación OFDM
        pulsos_modulados_ofdm = ifft(pulsos_ofdm); % Transformada inversa de Fourier

        % Agregar el prefijo cíclico
        pulsos_ofdm_cp = [pulsos_modulados_ofdm(end-cyclic_prefix_length+1:end, :); pulsos_modulados_ofdm]; % Agregar el prefijo cíclico

        % Agregar ruido a la señal recibida
        ofdm_awgn=awgn(pulsos_ofdm_cp,SNR(i),'measured');
       
        % Eliminar el prefijo cíclico
        senal_recibida_ofdm_sin_cp = ofdm_awgn(cyclic_prefix_length+1:end, :);

        % Demodulación OFDM
        pulsos_demodulados_ofdm = fft(senal_recibida_ofdm_sin_cp); % Transformada de Fourier
        
        bpsk_r=pskdemod(pulsos_demodulados_ofdm,Mod);  % Demodulación PSK
        bpsk_r_v = reshape(bpsk_r, size(x));            % convierto la matriz en vector columna
        isErr=x~=bpsk_r_v;                              % comparación be bits transmitidos con los recibidos
        BER(i,k)= nnz(isErr)/num_bits_s;
    end
end
Ber_promedio=mean(BER,2);
figure;
title("Rendimiento en canal awgn")
xlabel("SNR[dB]")
ylabel("Bit Error Rate")
semilogy (SNR, Ber_promedio)
grid on   