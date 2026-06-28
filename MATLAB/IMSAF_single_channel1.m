clear;
close all;

L = 1; % numero di altoparlanti (numero di canali)
M = 1; % numero di microfoni
I = 64; % numero di sottobande
D = 16; % fattore di decimazione
% V = 6*I+1; % lunghezza del filtro prototipo per il banco DFT
V = 129;
fs = 16000; % Frequenza di campionamento

K = 128; % lunghezza delle RIR (se la RIR vera è più lunga viene troncata)
offset = 100; % offset applicato alla RIR (utile se la RIR vera presenta molti 0 all'inizio)

P = 2;
% frameSize = 4096;
N = 300000; % lunghezza di x
mu_h = 0.5;
% delta_h = 0.1;
% delta_ap = 0.1;
delta_h = 1e-5; 
epsilon = 1e-5;
eps_reg = 1e-6;
delta_ap = 1e-4;

normalized_mis_flag = true; % Se impostato a true calcola il NM
norm_iteration_factor = 100;

fid = fopen("IR/IR1_1.f64");
H = fread(fid, 'double');
fclose(fid);
% Normalizzazione di H
H_original_peak = max(abs(H));
%H = H(offset+1:offset+K) / H_original_peak;
%H = H(offset+1:offset+K);
H = H/H_original_peak;

% RIR simulata per test
% H = zeros(8192, 1);
% H(10) = 1; 
% H(20) = -0.5; 
% H(40) = 0.2; 

K = length(H);
% if (D > 1)
%     Ki = ceil((K + V - 1) / D);
% else
%     Ki = ceil(K/D);
% end
Ki = ceil(K/D);
x = randn(N, 1);
% frameNumber = ceil(length(x)/frameSize);
s_subband = zeros(L*Ki, I); % matrice di I colonne, in cui la colonna i contiene
% le componenti di tutti i canali per la sottobanda i
s_ij = zeros(L*Ki, 1);
s_i_state = zeros(Ki, I);
y = zeros(K, L);
y_subband = zeros(floor(length(y)/D), I);
S_i_memory = zeros(L*Ki, P, I);
H_subband = zeros(I, Ki); % Dimensione I x Ki

% Idealmente i filtri devono tagliare a 1/I. Tuttavia, l'uso di filtri
% costruiti con fir1 limita l'esito della ricostruzione della RIR, portando
% il NM a non scendere sotto i -15dB

% Nota: usando un filtro equiripple il NM è più alto (circa -10dB) ma
% l'algoritmo converge dopo 1000 campioni

% prototype_dft_filter = fir1(V-1, 1/I);
 prototype_dft_filter = fir1(V-1, 1/I, kaiser(V, 12));
% prototype_dft_filter = firceqrip(V-1, 1/I, [0.05 0.03]);
% prototype_dft_filter = prototype_dft_filter / sum(prototype_dft_filter);

% Parametri rcosdesign
beta = 0.5;
span = 128;
%prototype_dft_filter = rcosdesign(beta, span, I, "sqrt");
%prototype_dft_filter = prototype_dft_filter / sqrt(I);
% Plot della risposta del filtro prototipo
figure;
freqz(prototype_dft_filter, 1, 512);

y = filter(H, 1, x);

%for frame=0:frameNumber-1
%input_frame = x(frame*frameSize+1:(frame+1)*frameSize);
%% Qui ci andrebbe il SFC q

s_subband = analysis_fb(x, prototype_dft_filter, I, D);
y_subband = analysis_fb(y, prototype_dft_filter, I, D);

K_len = size(s_subband, 1);
e_k = zeros(K_len, I);
u_ij = zeros(Ki,1);

norm_mis = zeros(1, floor(K_len/100));

% Calibrazione per il calcolo dei ritardi tramite impulso di test
impulse_len = K + 2 * V;
test_impulse = zeros(impulse_len, 1);
test_impulse(1) = 1;

% Si genera una delta di Dirac e la si fa passare attraverso i banchi di
% analisi e sintesi
sub_test = analysis_fb(test_impulse, prototype_dft_filter, I, D);
rec_test = synthesis_fb(sub_test, prototype_dft_filter, I, D);

% Salvando valore e posizione del picco si ottengono sia il ritardo
% introdotto dal banco filtri e sia il fattore di scala
[max_val, delay_calib] = max(real(rec_test));
delay_calib = delay_calib - 1; % Ritardo causato dal banco filtri
scale_factor = 1 / max_val; % Fattore di scala causato dal banco filtri

% Filtro di equalizzazione inverso usato per compensare le non idealità dei
% filtri

% Il filtro di equalizzazione viene ottenuto invertendo la trasformata
% dell'impulso fittizio transitato attraverso banco di analisi e sintesi.
% In ricostruzione, poi, si moltiplicherà la trasformata della RIR fullband
% e si ritornerà nel tempo

%FFT_len = 16384; % Potenza di 2 sufficientemente elevata
%G_fb = fft(rec_test, FFT_len);

% L'inversione della trasformata è regolarizzata da un epsilon per non
% farla divergere nei punti in cui G_fb -> 0. È importante calcolare la
% forma estesa e non 1 ./ (G_fb + eps_reg) per non alterare la fase.

%inv_G_fb = conj(G_fb) ./ (abs(G_fb).^2 + eps_reg);

% Array dei pesi sottobanda normalizzati (precalcolo dei pesi:
% implementazione statica).
sigma_i = sum(abs(s_subband).^2);
w = 1 ./ (sigma_i + epsilon);
w = w * I ./ sum(w);

for k=1:K_len
    for i=1:I
        s_i_state(:, i) = [s_subband(k, i); s_i_state(1:end-1, i)];
        s_ij = s_i_state(:, i);
        S_j = S_i_memory(:, :, i);
        % In MATLAB invece di inv() è consigliabile usare l'operatore \ che
        % risolve per x il sistema Ax=b, in quanto non è necessario
        % calcolare l'intera matrice inversa per poi moltiplicarla per un
        % vettore.

        % A rigore, si deve calcolare a = A^{-1}b con A = (S_j' * S_j) e b
        % = (S_j' * u_ij). Questo sistema può essere portato nella forma
        % Ax=b moltiplicando entrambi i membri per A, e ottenendo quindi
        % (S_j' * S_j) a = (S_j' * u_ij). Per funzionamento dell'operatore
        % \, quindi, a = (S_j' * S_j) \ (S_j'* u_ij).

        % NOTA: nell'articolo scrivono che a_i dipende da u ma poi lo
        % definiscono come "la proiezione ai minimi quadrati di s_i su
        % S_i", non è chiaro quale sia la formula corretta.

        % L'aggiunta di delta_ap*eye(P) somma una costante sulla diagonale
        % di S_j' * S_j, spostando i suoi autovalori dallo zero e evitando
        % instabilità dovuta al fatto che, per segnali correlati, S_j' *
        % S_j è quasi singolare
        a_i_col = (S_j' * S_j + delta_ap * eye(P)) \ (S_j' * s_ij);
        u_ij = s_ij - S_j * a_i_col;
        S_i_memory(:, :, i) = [s_ij, S_i_memory(:, 1:end-1, i)];
        e_k(k, i) = y_subband(k, i) - H_subband(i, :) * s_ij;
        norm_u = norm(u_ij)^2;

        % Nota: nel caso di x rumore bianco i pesi sono pressoché unitari
        H_subband(i, :) = H_subband(i, :) + mu_h * w(i) * e_k(k,i) * u_ij' / (norm_u + delta_h);
    end
    
    % Normalized misalignment calcolato ogni norm_iteration_factor campioni
    % per limitare la complessità
    if (normalized_mis_flag && mod(k, norm_iteration_factor) == 0)
        % Ricostruzione della H fullband:
        % per la ricostruzione corretta si fa passare l'impulso diviso in
        % sottobande attraverso le sottobande della RIR, e poi si
        % ricostruisce la RIR partendo da tali uscite

        % In alternativa, la H fullband può essere ricostruita facendo
        % upsampling della RIR stimata in sottobanda e calcolandone la 
        % convoluzione con i filtri di analisi e sintesi. 
        % Poiché la convoluzione è lineare, è possibile precalcolare la 
        % convoluzione dei filtri di analisi e sintesi e poi calcolare una 
        % sola convoluzione per ogni sottobanda.


        % Implementazione con filtraggio degli impulsi
        % y_sub_imp = zeros(size(sub_test));
        % for i_sub = 1:I
        %     y_sub_imp(:, i_sub) = filter(H_subband(i_sub, :), 1, sub_test(:, i_sub));
        % end
        % H_recon_raw = synthesis_fb(y_sub_imp, prototype_dft_filter, I, D);

        % Implementazione tramite convoluzione col filtro prototipo
        H_recon_raw = RIR_reconstruction(prototype_dft_filter, H_subband, I, D);

        % Equalizzazione della RIR
        % H_freq = fft(H_recon_raw, FFT_len);
        % H_recon_eq = real(ifft(H_freq .* inv_G_fb));
        % H_recon = H_recon_eq(1:K);
        H_recon = scale_factor * H_recon_raw(delay_calib + 1 : delay_calib + K);

        norm_mis(k/100) = 20*log10(norm(H-H_recon, 'fro')/norm(H, 'fro'));
    end
end

% Ricostruzione della H fullband
% y_sub_imp = zeros(size(sub_test));
% for i_sub = 1:I
%     y_sub_imp(:, i_sub) = filter(H_subband(i_sub, :), 1, sub_test(:, i_sub));
% end
% H_recon_raw = synthesis_fb(y_sub_imp, prototype_dft_filter, I, D);
H_recon_raw = RIR_reconstruction(prototype_dft_filter, H_subband, I, D);

% Compensazione dei ritardi e della scalatura introdotti dai filtri
% H_recon = scale_factor * H_recon_raw(delay_calib + 1 : delay_calib + K);
% Equalizzazione della RIR
% H_freq = fft(H_recon_raw, FFT_len);
% H_recon_eq = real(ifft(H_freq .* inv_G_fb));
% H_recon = H_recon_eq(1:K);
H_recon = scale_factor * H_recon_raw(delay_calib + 1 : delay_calib + K);

figure;
plot(H(100:228), 'LineWidth', 1.5);
hold on;
plot(H_recon(100:228), 'LineWidth', 1.5);
legend("Real RIR", "Estimated RIR");
title('RIR Time Domain Comparison');
xlabel('Samples');
ylabel('Normalized Amplitude');
grid on;

% Plot del NM
if (normalized_mis_flag)
    figure;
    plot((1:length(norm_mis))*D*(norm_iteration_factor/fs), norm_mis);
    xlabel('Iteration Time (s)');
    ylabel('NM (dB)');
    title('Normalized Misalignment');
    grid on;
end

% Plot degli spettri delle RIR
figure;
K_fft = length(H_recon);
H_true_fft = fft(H);
H_est_fft = fft(H_recon);
% Si mantengono solo le frequenze positive
half_len = floor(K_fft/2) + 1;
H_true_half = H_true_fft(1:half_len);
H_est_half = H_est_fft(1:half_len);
% Generazione del vettore delle frequenze fisiche in Hz da 0 a fs/2
freq_vec = (0 : half_len - 1) * (fs / K_fft);
mag_true_dB = 20 * log10(abs(H_true_half));
mag_est_dB = 20 * log10(abs(H_est_half));
plot(freq_vec, mag_true_dB, 'LineWidth', 1); 
hold on;
plot(freq_vec, mag_est_dB, 'LineWidth', 1);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('RIR Frequency Response Comparison');
legend("Real RIR", "Estimated RIR");
grid on;


%% Banchi di analisi e sintesi
function subband_signals = analysis_fb(x, hp, I, D)
    % x: segnale fullband di ingresso (vettore colonna)
    % hp: filtro prototipo di lunghezza V
    % I: numero di sottobande
    % D: fattore di decimazione
    
    N = length(x);
    
    % Numero di campioni nel dominio subband decimato
    K_len = floor((N - 1) / D) + 1;
    subband_signals = zeros(K_len, I);
    n_idx = (0:N-1)';
    for i = 0:I-1
        % Modulazione del segnale di ingresso
        x_mod = x .* exp(-1i * 2 * pi * i * n_idx / I);
        % Filtraggio Passa-Basso
        v_i = filter(hp, 1, x_mod);
        % Decimazione
        subband_signals(:, i+1) = v_i(1:D:end);
    end
end

function x_recon = synthesis_fb(subband_signals, hp, I, D)
    % subband_signals: matrice di dimensione (K_len x I)
    % hp: filtro prototipo di lunghezza V
    % I: numero di sottobande
    % D: fattore di decimazione
    
    K_len = size(subband_signals, 1);
    N_recon = (K_len - 1) * D + 1;
    x_recon = zeros(N_recon, 1);
    n_idx = (0:N_recon-1)';
    for i = 0:I-1
        % Upsampling
        w_i = zeros(N_recon, 1);
        w_i(1:D:end) = subband_signals(:, i+1);
        % Filtraggio Passa-Basso (interpolazione)
        z_i = filter(hp, 1, w_i);
        % Modulazione inversa
        s_i = z_i .* exp(1i * 2 * pi * i * n_idx / I);
        
        x_recon = x_recon + real(s_i);
    end
    % Compensazione del guadagno
    % x_recon = x_recon / D;
end

%% Ricostruzione della RIR fullband
function w_full = RIR_reconstruction(hp, H_subband, I, D)
    % Funzione che ricostruisce la RIR fullband partendo dai valori in
    % sottobande e sfruttando il filtro prototipo.

    % hp: filtro prototipo di lunghezza N_filt
    % H_subband: matrice IxKi che contiene le componenti subband della RIR
    % I: numero di sottobande
    % D: fattore di decimazione
    % delay: delay introdotto dai banchi filtro

    % Assicura che il hp sia un vettore colonna
    hp = hp(:); 
    N_filt = length(hp);
    Ki = size(H_subband, 2);
    
    L_upsampled = (Ki - 1) * D + 1;
    L_total = N_filt + L_upsampled + N_filt - 2;
    
    w_full = zeros(L_total, 1);
    t_total_idx = (0:L_total-1)';
    
    % Si può precalcolare la convoluzione del filtro prototipo con sé
    % stesso e poi fare la convoluzione di questo segnale per le componenti
    % in sottobanda
    combined_filter = conv(hp, hp);
    for i = 0:I-1
        % Upsampling
        w_i = H_subband(i+1, :).';
        w_i_upsampled = zeros(L_upsampled, 1);
        w_i_upsampled(1:D:end) = w_i;
        
        % Doppia convoluzione reale in banda base (che corrisponde
        % all'applicazione del filtro di analisi e sintesi)
        % w_base = conv(conv(hp, w_i_upsampled), hp);
        w_base = conv(w_i_upsampled, combined_filter);
        
        % Modulazione complessa del risultato finale fullband
        t_i = w_base .* exp(1i * 2 * pi * i * t_total_idx / I);
        
        % Accumulo
        w_full = w_full + real(t_i);
    end
    
    % Compensazione del fattore di scala
    w_full = w_full / D;
end

