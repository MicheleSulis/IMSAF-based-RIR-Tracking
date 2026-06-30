function plot_RIR(RIR, K, fs, offset)
figure;
H_fft = fft(RIR(offset+1:offset+K)/max(abs(RIR(offset+1:offset+K))));
half_len = floor(K/2) + 1;
H_half = H_fft(1:half_len);
% Generazione del vettore delle frequenze fisiche in Hz da 0 a fs/2
freq_vec = (0 : half_len - 1) * (fs / K);
mag_est_dB = 20 * log10(abs(H_half));
plot(freq_vec, mag_est_dB, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('RIR Frequency Response');
