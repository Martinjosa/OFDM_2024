clear
close all
clc

% Fijar la semilla del generador de números aleatorios
seed = 4;  % Se puede elegir cualquier número
rng(seed);

%======================== Inicio de variables ============================%
num_bits = 100000;              % Número de bits a generar
num_subportadoras = 64;         % Número de subportadoras en OFDM
cyclic_prefix_length = 16;      % Longitud del prefijo cíclico en OFDM
num_pilotos = 8;                % Número de pilotos
M = 2;                          % Bits por símbolo [1 2 4]
k = log2(M);                    % Bits por símbolo
delta_f = 1;                     % Separación en frecuencia entre subportadoras
delta_f_c = 20;                 % Ancho de banda de coherencia según entorno

% Ajustar el número de bits para que sea un múltiplo del número de subportadoras
num_bits_s = ceil(num_bits / (num_subportadoras - num_pilotos)) * (num_subportadoras - num_pilotos);

%======================== TRANSMISOR OFDM ========================%
%=================================================================%

pulsos_digitales = randi([0, 1], num_bits_s, 1);               % Genero el vector columna, bits en paralelo
bpsk_signal = pskmod(pulsos_digitales, M);                     % Modulacion BPSK

% Generar símbolos OFDM serial to paralelo
num_symbols = ceil(num_bits_s / (num_subportadoras - num_pilotos)); % Número de simbolos OFDM necesarios
pulsos_ofdm = reshape(bpsk_signal, num_subportadoras - num_pilotos, num_symbols); % Agrupar los bits en símbolos OFDM

% Insertar pilotos
indices_pilotos = round(linspace(1, num_subportadoras, num_pilotos)); % Índices de las subportadoras piloto de forma aleatoria
%indices_pilotos = 1:num_pilotos:num_subportadoras; % Índices de las subportadoras piloto de forma estática
secuencias_piloto = ones(length(indices_pilotos), num_symbols); % Secuencias de pilotos (puede ser cualquier secuencia conocida)

pulsos_ofdm_con_pilotos = zeros(num_subportadoras, num_symbols);
pulsos_ofdm_con_pilotos(indices_pilotos, :) = secuencias_piloto; % Insertar pilotos en los índices correspondientes
pulsos_ofdm_con_pilotos(~ismember(1:num_subportadoras, indices_pilotos), :) = pulsos_ofdm; % Insertar datos en las demás subportadoras

% Modulación OFDM
pulsos_modulados_ofdm = ifft(pulsos_ofdm_con_pilotos, num_subportadoras); % Transformada inversa de Fourier

% Agregar el prefijo cíclico
pulsos_ofdm_cp = [pulsos_modulados_ofdm(end-cyclic_prefix_length+1:end, :); pulsos_modulados_ofdm]; % Agregar el prefijo cíclico

% Convertir a serie para el canal
pulsos_modulados_ofdm_serie = pulsos_ofdm_cp(:);

% FIN TRANSMISOR OFDM

        %================= CANAL RAYLEIGH + AWGN =========================%
        %=================================================================%

% Parámetros del canal de comunicación
EbNo = 25;                      % Energia de bits a ruido
SNR_dB = EbNo + 10*log10(k);    % Relación señal-ruido en dB
SNR = 10^(SNR_dB / 10);         % Relación señal-ruido en escala lineal
ruido = sqrt(1 / (2*SNR));      % Ruido AWGN

% Simular el canal de Rayleigh
rayChan = comm.RayleighChannel( ...
    'SampleRate', 1, ...
    'PathDelays', [0 1.5e-5 3.2e-5], ...
    'AveragePathGains', [0 -2 -10], ...
    'NormalizePathGains', true, ...
    'MaximumDopplerShift', 0);

% Pasar la señal OFDM a través del canal
OFDM_Ray = rayChan(pulsos_modulados_ofdm_serie);

% Añadir ruido AWGN
ofdm_awgn = OFDM_Ray + ruido * (randn(size(OFDM_Ray)) + 1i * randn(size(OFDM_Ray)));
%ofdm_awgn = awgn(OFDM_Ray,SNR_dB,ruido);


 %==================== RECEPTOR OFDM ==============================%
 %=================================================================%

% Convertir de nuevo a forma paralela
ofdm_awgn_parallel = reshape(ofdm_awgn, num_subportadoras + cyclic_prefix_length, num_symbols);

% Eliminar el prefijo cíclico
senal_recibida_ofdm_sin_cp = ofdm_awgn_parallel(cyclic_prefix_length+1:end, :);

% Demodulación OFDM
pulsos_demodulados_ofdm = fft(senal_recibida_ofdm_sin_cp, num_subportadoras); % Transformada de Fourier

% Estimación del canal usando los pilotos (MMSE)
% Calcular la matriz de autocorrelación del canal y la matriz de ruido
rho = exp(-delta_f/ delta_f_c); % Referído a la correlación espacial
Rhh = rho .^ abs(repmat((1:num_pilotos)', 1, num_pilotos) - repmat(1:num_pilotos, num_pilotos, 1));  % Autocorrelación del canal  
Rnn = (1 / SNR) * eye(num_pilotos); % Autocorrelación del ruido
%Rhh = toeplitz(rho.^(0:8-1));
% Estimación del canal en los índices de los pilotos
H_est_pilotos = pulsos_demodulados_ofdm(indices_pilotos, :) ./ secuencias_piloto; 

% Estimación MMSE del canal
H_mmse = (Rhh / (Rhh + Rnn)) * H_est_pilotos;

% Interpolación lineal para obtener la estimación del canal completa
H_interpolado = interp1(indices_pilotos, H_mmse, 1:num_subportadoras, 'linear', 'extrap');

% Ecualización
yEq = pulsos_demodulados_ofdm ./ H_interpolado;

% Extraer datos de las subportadoras que no son pilotos
datos_recibidos = yEq(~ismember(1:num_subportadoras, indices_pilotos), :);

% Convertir de nuevo a serie
pulsos_demodulados_ofdm_serie = reshape(datos_recibidos, 1, []);

bpsk_r = pskdemod(pulsos_demodulados_ofdm_serie, M);  % Demodulación PSK

% FIN RECEPTOR OFDM

pulsos_digitales_fila = pulsos_digitales.'; % paso a vector fila los bits transmitidos

% Visualizar los primeros 15 símbolos transmitidos
disp('Símbolos transmitidos:');
disp(pulsos_digitales_fila(1:15));

% Visualizar los primeros 15 símbolos demodulados
disp('Símbolos demodulados:');
disp(bpsk_r(1:15));

% Cálculo de la tasa de error de bit (BER)
ber_otra = sum(pulsos_digitales_fila ~= bpsk_r) / num_bits_s;
disp(['Tasa de error de bit (BER): ', num2str(ber_otra)]);
[numErrors, ber] = biterr(pulsos_digitales_fila, bpsk_r);
fprintf('\nLa tasa de error de bits de codificación binaria es %5.2e, basada en %d errores.\n', ber, numErrors);

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
ylim([-2.5, 2.5]); % Establecer los límites del eje y
