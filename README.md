# OFDM_2024
Este repositorio es para documentar los avances de la simulación de una sistema OFDM. El objetivo es constuir un transmisor OFDM, esta señal modulada pasarla por un canal con ruido gauseano AWGN y desvanecimiento por trayectos multiples Rayleigh y en el extremo reseptor hacer la demodulacion OFDM y la correcta ecualización para contrarestar los efectos del canal. Asi mismo, al sistema se lo someterá a la simulacion Montecarlo para determinar el rendimiento y la BER con diferentes niveles en la relacion señal-ruido SNR.

La idea es ir contruyendo el sistema de apoco:
  1- Se crea un sistema BPSK.
  2- se le agrega ruido AWGN.
  3- Simulacion Montecarlo.
  4- Sobre el sistema BPSK se agrega OFDM con prefijo ciclico.
  5- Se agrega el canal Rayleigh basico y se agrega un ecualizador basico.
  6- Simulación Montecarlo.
  7- Se crea el canal Rayleigh como objeto, se introducen portadoras pilotoas en el transmisor y se ecualiza MSE en el recepto.
  8- Simulación Montecarlo.
  
