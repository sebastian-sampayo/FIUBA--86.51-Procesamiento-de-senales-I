% Test de WLS
% versión 1

clear all;
% close all;

M = 101;
A = [ones(40,1); zeros(60,1)];
W = ones(100,1);
flg_type = 4;
% A_bands = [0 0 1 1 0 0 .5 .5 0 0];
% f_bands = [0 .1 .2 .3 .4 .55 .6 .8 .85 1];
A_bands = [1 1];
f_bands = [.03 1];
delta0 = .01;
delta1 = .1;
delta2 = .01;
delta3 = .05;
delta4 = .01;
% W_bands = [1 1 1 1];
% W_bands = [1/delta0 1/delta1 1/delta2 1/delta3 1/delta4];
W_bands = [1];

% h_v1 = WLS_v1(M, A, W, flg_type);
h_v2 = WLS(M, f_bands, A_bands, W_bands, flg_type);
h_firls = firls(M, f_bands, A_bands, W_bands, 'hilbert');
h_pm = firpm(M, f_bands, A_bands, W_bands, 'hilbert');
% Ventaneo del h ideal
n = (0:M)';
h_ideal = zeros(M, 1);
h_ideal = 1./(pi*(n-M/2)) .* (1- cos((n-M/2)*pi)) .*(n~=M/2);
h_ventaneo = h_ideal .* hamming(M+1);

% wvtool(h_v2, h_firls);
fvtool(h_v2,1, h_firls,1, h_ventaneo,1, h_pm,1,...
    'MagnitudeDisplay','magnitude');
% legend('WLS', 'firls');
% 
% nfft = 2^10;
% H_v2 = fft(h_firls.*hamming(M+1)', nfft);
% arg_v2 = angle(H_v2);
% figure
% plot(arg_v2)
% figure
% hold all
% stem(-h_v2)
% stem(h_firls)

