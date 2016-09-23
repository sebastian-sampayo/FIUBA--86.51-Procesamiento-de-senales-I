% -------------------------------------------------------------------------- %
% Facultad de Ingeniería de la Universidad de Buenos Aires
% Procesamiento de Señales I
% Trabajo Práctico 1: Diseño de filtro con fase lineal
%   - Transformador de Hilbert para codificación de voz -
% 1° Cuatrimestre de 2015
%
% Sampayo, Sebastián Lucas
% Padrón: 93793
% e-mail: sebisampayo@gmail.com
%
% Script principal - MATLAB
% -------------------------------------------------------------------------- %

close all;
clear all;

% Nombre del archivo de entrada (Audio codificado)
input_file_name = 'pe_revuelto.wav';
output_file_name = 'pe_decoded.wav';

% Frecuencia de referencia para el "espejado" $\Omega_r$
Omega_r = 4100; % [Hz]
% Frecuencia de paso, a partir de la cual comienza 
% la banda de paso para el Transformador de Hilbert
Omega_p = 100; % [Hz]
% Para este caso, se eligieron 100Hz, considerando el ligero aliasing que
% ocurre en este caso en bajas frecuencias. Por otro lado, la mayor parte de 
% la información del habla se encuentra arriba de esta frecuencia.


% Cargo el archivo de entrada
[x, fs] = audioread(input_file_name);
x = x(:,1); % Me quedo con un solo canal (en caso de que sea estéreo)
Ts = 1/fs;
% Frecuencia de referencia para el "espejado" en dominio discreto $\omega_r$
wr = Omega_r * 2*pi/fs;
% Frecuencia de paso, a partir de la cual comienza 
% la banda de paso para el Transformador de Hilbert en el dominio discreto
wp = Omega_p * 2*pi/fs; 

% --------------------------------------------------------------------------- %
% ---- Diseño del Transformador de Hilbert FIR FLG ----
% --------------------------------------------------------------------------- %
M = 187;
A = [1 1];
f = [wp/pi 1]; % Divido por \pi para normalizar
delta = .01;
W = [1/delta]; 
flg_type = 4;
%   Se desea que el filtro sea antisimétrico y valga 'cero' en omega=0 y 'uno' en \pi.
% Los filtros tipo I y II son simétricos; y el tipo III fuerza un cero en pi, por lo
% tanto estas configuraciones fueron descartadas. Se optó entonces por un Tipo IV que
% a su vez, asegura un 'cero' en omega=0 y antisimetría.
%   El delta es una formalidad, ya que al poner un único valor constante
% para la función peso (W), el resultado es igual a poner cualquier constante.
% Por lo tanto, la especificación de error en la banda de paso deberá
% vigilarse con la simulación, modificando el orden M hasta alcanzar el error
% deseado con algún filtro.

% Por Cuadrados mínimos
h_wls = WLS(M, f, A, W, flg_type);
% Por Cuadrados mínimos con la función 'firls' de MATLAB
h_firls = firls(M, f, A, W, 'hilbert')';
% Por algoritmo Parks-McClellan
h_pm = firpm(M, f, A, W, 'hilbert')';
% Por ventaneo de la respuesta al impulso del filtro ideal
n = (0:M)';
h_ideal = zeros(M, 1);
h_ideal = 1./(pi*(n-M/2)) .* (1- cos((n-M/2)*pi)) .*(n~=M/2);
window = hanning(M+1);
h_ventaneo = h_ideal .* window;

% wvtool(h_v2, h_firls);
fvtool(h_wls,1, h_pm,1, h_ventaneo,1, h_firls,1,...
    'MagnitudeDisplay','magnitude');

% Se alcanza la especificación utilizando 
% Omega_p = 100; % [Hz]
% M = 187;
% A = [1 1];
% f = [wp/pi 1];
% W = [1/delta]; 

% Y relajando más la banda de transición se obtienen ordenes menores:
% Omega_p = 200; % [Hz]
% M = 97;

% Omega_p = 300; % [Hz]
% M = 63;

%   También se puede destacar que el filtro de Parks_McClellan alcanza la 
% especificación con un orden ligeramente menor que el filtro de WLS.
% No obstante, la energía del error para el primero se mantiene constante
% a lo largo de la banda de paso, mientras que para el segundo, el error
% se extingue muy rápidamente, alcanzando mayor precisión en la mayor parte
% de la banda. Por esta razón, se eligió el segundo para la realización
% del sistema de codificación.
%   Adicionalmente, se incluye una solución por ventaneo de Hanning, en la cual
% se observa que es el diseño cuya energía de error a lo largo de la banda
% pareciera ser menor que en ambos casos anteriores (no calculado, pero se 
% desprende del gráfico).

% Filtro elegido h:
h = h_wls;
% --------------------------------------------------------------------------- %

% --------------------------------------------------------------------------- %
% ---- Sistema de codificación ----
% --------------------------------------------------------------------------- %
% Filtro retardo
h_delay = sinc(n-M/2);

% Le aplico los filtros obtenidos a la entrada 'x'
x_hilbert = conv(x, h);
x_delay = conv(x, h_delay);

N_x = length(x);
n = (0:(M+N_x-1))';

% Multiplicación por seno y coseno (modulación)
z1 = x_hilbert .* sin(wr*n);
z2 = x_delay .* cos(wr*n);

% Combinación final
y = -z1 - z2;
% --------------------------------------------------------------------------- %

% --------------------------------------------------------------------------- %
% ---- Gráficos ----
% --------------------------------------------------------------------------- %
% Espectros
nfft = 2^11;
X_fft = fft(x.*hamming(N_x), nfft);
X_fft = X_fft(1:nfft/2);

N_y = length(y);
Y_fft = fft(y.*hamming(N_y), nfft);
Y_fft = Y_fft(1:nfft/2);

% window = hamming(M+1);
window = ones(M+1,1);

H_fft = fft(h.*window, nfft);
H_fft = H_fft(1:nfft/2);

H_wls_fft = fft(h_wls.*window, nfft);
H_wls_fft = H_wls_fft(1:nfft/2);

H_pm_fft = fft(h_pm.*window, nfft);
H_pm_fft = H_pm_fft(1:nfft/2);

H_ventaneo_fft = fft(h_ventaneo.*window, nfft);
H_ventaneo_fft = H_ventaneo_fft(1:nfft/2);

H_firls_fft = fft(h_firls.*window, nfft);
H_firls_fft = H_firls_fft(1:nfft/2);

w_n = linspace(0,1, nfft/2); % frecuencias normalizadas en \pi

figure(2)
title('Respuesta en frecuencia de señales de entrada y salida');
hold all;
plot(w_n, abs(X_fft))
plot(w_n, abs(Y_fft))
legend('Espectro de la señal de entrada (codificada)',...
        'Espectro de la señal de salida (decodificada)');
grid on;
xlabel('Frecuencia (Normalizada en \pi) [rad/\pi]');
% print('-dpng', 'Espectros_señales_in_out.png');

figure(3)
title('Respuesta en frecuencia de filtros de Hilbert diseñados');
hold all;
plot(w_n, abs(H_wls_fft))
plot(w_n, abs(H_pm_fft))
legend('Cuadrados Mínimos Ponderados',...
        'Algoritmo Parks-McClellan');
plot([0 1], [1+delta 1+delta], 'k--')
plot([0 1], [1-delta 1-delta], 'k--')
plot([wp wp]/pi, [0 1.1], 'k--')
text(-.05, 1+delta, '1+\delta')
text(-.05, 1-delta, '1-\delta')
text(wp/pi, -.05, '\omega_p')
grid on;
xlabel('Frecuencia (Normalizada en \pi) [rad/\pi]');
ylim([0 1.1]);
% print('-dpng', 'Filtros_de_Hilbert.png');

figure(4)
title('Comparación de la función "WLS" y "firls"')
hold all;
plot(w_n, abs(H_wls_fft))
plot(w_n, abs(H_firls_fft))
legend('WLS', 'firls');
grid on;
xlabel('Frecuencia (Normalizada en \pi) [rad/\pi]');
ylim([0 1.1]);
str = sprintf('max(|h_{wls}(n) - h_{firls}(n)|) = %f', max(abs(h_wls-h_firls)));
text(.1,.8,str);
% print('-dpng', 'WLS_vs_firls.png');

figure(5)
title('Respuesta al impulso del filtro de Hilbert diseñado')
hold all;
n = (0:M)';
stem(n,h)
grid on;
text(M/2, -.01, 'M/2');
% print('-dpng', 'Respuesta_al_impulso_Hilbert.png');


% --------------------------------------------------------------------------- %
% ---- Salida ----
% --------------------------------------------------------------------------- %
% Normalizo a picos del 80% del volumen máximo posible
y = .8*y/max(abs(y));

% Imprimo salida en archivo de audio
% audiowrite(output_file_name, y, fs);
