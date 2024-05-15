%========================= Borrar datos ==================================%
%=========================================================================%
% Este ejemplo de canal se saco de https://www.youtube.com/watch?v=NCd1Cs2WMf0&t=1303s

clear
close all
clc

%======================== Inicio de variables ============================%
num_bits = 100000;                                  % 100 mil bits
SNR=[1:30];
num_subportadoras = 64; % Número de subportadoras en OFDM
cyclic_prefix_length = 16; % Longitud del prefijo cíclico en OFDM
Mod = 2; % Bits por símbolo [1 2 4]

% Ajustar el número de bits para que sea un múltiplo del número de subportadoras
num_bits_s = ceil(num_bits / num_subportadoras) * num_subportadoras;

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

        h = [0.8; zeros(7,1); -0.5; zeros(7,1); 0.34]; %Canal Rayleigh
        OFDM_Ray = filter(h,1,pulsos_ofdm_cp); %Pasa la señal por el canal

        
        % Agregar ruido a la señal recibida
        ofdm_awgn=awgn(OFDM_Ray,SNR(i),'measured');
       
        % Eliminar el prefijo cíclico
        senal_recibida_ofdm_sin_cp = ofdm_awgn(cyclic_prefix_length+1:end, :);

        % Demodulación OFDM
        pulsos_demodulados_ofdm = fft(senal_recibida_ofdm_sin_cp); % Transformada de Fourier
        
        % Respuesta en frecuencia del canal estimada
        H = fftshift(fft(h, num_subportadoras));

        %Ecualizacion
        yEq = pulsos_demodulados_ofdm ./ H;

        
        bpsk_r=pskdemod(yEq,Mod);  % Demodulación PSK
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