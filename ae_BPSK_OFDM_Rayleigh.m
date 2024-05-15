%========================= Borrar datos ==================================%
%=========================================================================%
% Este ejemplo de canal se saco de https://www.youtube.com/watch?v=NCd1Cs2WMf0&t=1303s

clear
close all
clc

% Parámetros de la fuente de pulsos digitales
num_bits = 100000;              % Número de bits a generar
num_subportadoras = 64;         % Número de subportadoras en OFDM
cyclic_prefix_length = 16;      % Longitud del prefijo cíclico en OFDM
M = 2;                          % Bits por símbolo [1 2 4]
k = log2(M);                    %bits por simbolo

% Ajustar el número de bits para que sea un múltiplo del número de subportadoras
num_bits_s = ceil(num_bits / num_subportadoras) * num_subportadoras;


 pulsos_digitales = randi([0,1],num_bits_s,1);                % Genero el vector columna, bits en paralelo
   
 bpsk_signal = pskmod(pulsos_digitales,M);                   % Modulacion bpsk
        
% Generar símbolos OFDM serial to paralelo
num_symbols = ceil(num_bits_s / num_subportadoras); % Número de simbolos OFDM necesarios
pulsos_ofdm = reshape(bpsk_signal, num_subportadoras, num_symbols); % Agrupar los bits en símbolos OFDM (las 64 portadoras en paralelo)
        
% Modulación OFDM
pulsos_modulados_ofdm = ifft(pulsos_ofdm); % Transformada inversa de Fourier

% Agregar el prefijo cíclico
pulsos_ofdm_cp = [pulsos_modulados_ofdm(end-cyclic_prefix_length+1:end, :); pulsos_modulados_ofdm]; % Agregar el prefijo cíclico


% Parámetros del canal de comunicación
EbNo = 35;                      % Energia de bits a ruido
SNR_dB = EbNo + 10*log10(k);    % Relación señal-ruido en dB
SNR = 10^(SNR_dB / 10);         % Relación señal-ruido en escala lineal
ruido = sqrt(1 / SNR);     

h = [0.8; zeros(7,1); -0.5; zeros(7,1); 0.34];
OFDM_Ray = filter(h,1,pulsos_ofdm_cp);

   
ofdm_awgn = awgn(OFDM_Ray,SNR_dB,ruido);
   

% Eliminar el prefijo cíclico
senal_recibida_ofdm_sin_cp = ofdm_awgn(cyclic_prefix_length+1:end, :);

% Demodulación OFDM
pulsos_demodulados_ofdm = fft(senal_recibida_ofdm_sin_cp); % Transformada de Fourier

% Respuesta en frecuencia del canal estimada
H = fftshift(fft(h, num_subportadoras));

%Ecualizacion
yEq = pulsos_demodulados_ofdm ./ H;


 pulsos_demodulados_ofdm_serie = reshape(yEq, 1, []);
 
 bpsk_r=pskdemod(pulsos_demodulados_ofdm_serie,M);  % Demodulación PSK
 
 pulsos_digitales_fila = pulsos_digitales.'; % paso a vector fila
 bpsk_r_vector = reshape(bpsk_r, 1, []);% paso a vector fila
 
  % Visualizar los primeros 10 símbolos transmitidos
  disp('Símbolos transmitidos:');
  disp(pulsos_digitales_fila(1:15));
  
  % Visualizar los primeros 10 símbolos demodulados
  disp('Símbolos demodulados:');
  disp(bpsk_r(1:15));
 
  % Cálculo de la tasa de error de bit (BER)
  
  ber_otra = sum(pulsos_digitales_fila ~= bpsk_r) / num_bits_s;
   disp(['Tasa de error de bit (BER): ', num2str(ber_otra)]);
   [numErrors,ber] = biterr(pulsos_digitales_fila,bpsk_r);
   fprintf( '\nLa tasa de error de bits de codificación binaria es %5.2e, basada en %d errores.\n' , ... 
       ber,numErrors)
 
 % Crear un diagrama de constelación de los símbolos demodulados
 figure;
 
 scatter(real(pulsos_demodulados_ofdm_serie), imag(pulsos_demodulados_ofdm_serie), '.');
 title('Diagrama de constelación de los símbolos demodulados');
 xlabel('Parte Real');
 ylabel('Parte Imaginaria');
 axis equal; % Para mantener la escala de los ejes iguales
 grid on; % Mostrar una cuadrícula en el gráfico
 
% Ampliar el rango del eje x
 xlim([-2.5, 2.5]); % Establecer los límites del eje x
