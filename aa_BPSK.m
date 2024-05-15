%========================= Borrar datos ==================================%
%=========================================================================%
close all
clear

% Parámetros de la fuente de pulsos digitales
num_bits = 100000; % Número de bits a generar

% Generar pulsos digitales aleatorios
pulsos_digitales = randsample([0, 1], num_bits, true); % Genera 0 y 1

% Parámetros de modulación PSK
fase_0 = 0; % Fase del símbolo 0
fase_1 = pi; % Fase del símbolo 1
%frecuencia_portadora = 1e9; % Frecuencia de la portadora (1 GHz)
%tiempo_bit = 1 / frecuencia_portadora; % Duración de cada bit en segundos

% Modulación PSK
pulsos_modulados = exp(1i * (pulsos_digitales * fase_1 + (1 - pulsos_digitales) * fase_0)); % Modulación PSK


% Parámetros del canal de comunicación
SNR_dB = 15; % Relación señal-ruido en dB
SNR = 10^(SNR_dB / 10); % Relación señal-ruido en escala lineal
ruido = sqrt(1 / SNR);

senal_recibida = awgn(pulsos_modulados,SNR_dB,ruido);  % Generar ruido gaussiano
% Simular la recepción de la señal a través de la antena receptora
%senal_recibida = pulsos_modulados + ruido;


% Demodulación PSK
  pulsos_demodulados = zeros(size(senal_recibida));
 %pulsos_demodulados(angle(senal_recibida) > fase_1/2) = 1;
% Calcular el ángulo de la señal recibida
angulos = atan2(imag(senal_recibida), real(senal_recibida));

% Decidir si cada símbolo demodulado es 0 o 1 en función del ángulo
pulsos_demodulados(angulos > pi/2 | -pi/2 > angulos) = 1;


% Visualizar los primeros 20 símbolos transmitidos
disp('Símbolos transmitidos:');
disp(pulsos_digitales(1:20));

% Visualizar los primeros 20 símbolos demodulados
disp('Símbolos demodulados:');
disp(pulsos_demodulados(1:20));

% Cálculo de la tasa de error de bit (BER)
ber_otra = sum(pulsos_digitales ~= pulsos_demodulados) / num_bits;
disp(['Tasa de error de bit (BER): ', num2str(ber_otra)]);
[numErrors,ber] = biterr(pulsos_digitales,pulsos_demodulados);
fprintf( '\nLa tasa de error de bits de codificación binaria es %5.2e, basada en %d errores.\n' , ... 
    ber,numErrors)

% Crear un diagrama de constelación de los símbolos demodulados
figure;
scatter(real(senal_recibida), imag(senal_recibida), '.');
title('Diagrama de constelación de los símbolos demodulados');
xlabel('Parte Real');
ylabel('Parte Imaginaria');
axis equal; % Para mantener la escala de los ejes iguales
grid on; % Mostrar una cuadrícula en el gráfico

% Ampliar el rango del eje x
xlim([-2.5, 2.5]); % Establecer los límites del eje x
