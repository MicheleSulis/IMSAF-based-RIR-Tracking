clear;
close all;

L = 2; % numero di altoparlanti (numero di canali)
M = 2; % numero di microfoni
M_plotted = 1;
L_plotted = 2;
I = 8; % numero di sottobande
D = 2; % fattore di decimazione
P = 2; % P = 0: nessuna decorrelazione
% N = 200000; % lunghezza di x
mu_h = 0.32;
delta_h = 1e-3;  
delta_ap = 1e-1;
K = 128; % lunghezza delle RIR (se la RIR vera è più lunga viene troncata)
offset = 100; % offset applicato alla RIR (utile se la RIR vera presenta molti 0 all'inizio)
fs = 16000; % Frequenza di campionamento (necessaria per la modulazione di fase)
N = fs * 30; % Permette di scegliere N in secondi a partire da fs

% Flag di test:
% - 'true': genera segnali di test completamente incorrelati su tutti i
% canali. Utile per mostrare che la struttura dell'algoritmo di tracking è
% corretta
% - 'false': viene usato un segnale reale, inviato sui vari canali.
x_test = false;

% Flag per scegliere il segnale di ingresso quando x_test è false
x_type = 2;
% 1: segnale mono replicato su tutti i canali, caricato dalla libreria di
% MATLAB
% 2: segnale stereo inviato a canali separati, preso da EBU SQUAM

normalized_mis_flag = true; % Se impostato a true calcola il NM
norm_iteration_factor = 100;

% V = 128;
% prototype_dft_filter = fir1(V-1, 1/I);
% prototype_dft_filter = fir1(V-1, 1/I, kaiser(V, 12));
% prototype_dft_filter = firceqrip(V-1, 1/I, [0.05 0.03]);
% prototype_dft_filter = prototype_dft_filter / sum(prototype_dft_filter);

% Parametri rcosdesign
beta = 0.6;
span = 16;
prototype_dft_filter = rcosdesign(beta, span, I, "sqrt");
prototype_dft_filter = prototype_dft_filter / sqrt(I);
V = length(prototype_dft_filter);

% Plot della risposta del filtro prototipo
figure;
freqz(prototype_dft_filter, 1, 512);

H = zeros(K, M, L);
file_index = 1;
for m=1:M
    for l=1:L
        fid = fopen("IR\IR1_"+string(file_index)+".f64");
        RIR = fread(fid, 'double');
        fclose(fid);
        H_original_peak = max(abs(RIR(offset+1:offset+K)));
        H(:, m, l) = RIR(offset+1:offset+K) / H_original_peak;
        % temp_H = zeros(K, 1);
        % temp_H(20) = 1;
        % H(:, m, l) = temp_H;
        file_index = file_index + 1;
    end
end
Ki = ceil((K + V)/D);
K_len = floor((N - 1) / D) + 1;
N_recon = (K_len - 1) * D + 1;
H_hat = zeros(M, L*Ki);
s_subband_decorr = zeros(K_len, I, L);
s_subband_base = zeros(K_len, I, L);

% Generazione del segnale di eccitazione (a singolo canale, viene reso
% multicanale dal SFC e poi decorrelato)
% x_mono = randn(N, 1);

if (x_type == 1)
    load handel.mat; x_audio = y; fs_audio = Fs;
    if fs_audio ~= fs
        x_audio = resample(x_audio, fs, fs_audio);
    end
    if length(x_audio) >= N
        x_mono = x_audio(1:N);
    else
        x_mono = repmat(x_audio, ceil(N/length(x_audio)), 1);
        x_mono = x_mono(1:N);
    end
    x_mono = x_mono / std(x_mono); 
    x_mono_subband = analysis_fb(x_mono, prototype_dft_filter, I, D);
    s_subband_base = repmat(x_mono_subband, 1, 1, L);
elseif (x_type == 2)
    [x_audio, fs_audio] = audioread("SampleAudio\69.flac");
    
    if fs_audio ~= fs
        x_audio = resample(x_audio, fs, fs_audio);
    end
    
    if length(x_audio) >= N
        x_stereo = x_audio(1:N, :);
    else
        x_stereo = repmat(x_audio, ceil(N/length(x_audio)), 1);
        x_stereo = x_stereo(1:N, :);
    end
    
    % Normalizzazione della potenza
    x_stereo(:, 1) = x_stereo(:, 1) / std(x_stereo(:, 1));
    x_stereo(:, 2) = x_stereo(:, 2) / std(x_stereo(:, 2));
    
    % Segnali stereo
    for l = 1:L
        s_subband_base(:,:,l) = analysis_fb(x_stereo(:, l), prototype_dft_filter, I, D);
    end
end

% Rumore marrone: sono riportate le prestazioni in termini di NM
% Canali scorrelati
%   -55dB dopo 383 iterazioni @ P=0;
%   -55dB dopo  82 iterazioni @ P=2.
% Canali correlati
%    -8dB dopo 470 iterazioni @ P=0;
%    -8dB dopo 326 iterazioni @ P=2.
% In entrambi i casi i floor sono paragonabili tra di loro e rispetto a
% x_mono = rumore bianco. Le convergenze sono rispettivamente 4.67 e 1.44
% volte più veloci passando da P=0 a P=2, provando, come affermato 
% nell'articolo che lo sbiancamento migliora la velocità di convergenza e 
% non il floor del NM.
x_mono = filter(1,[1 -0.9],x_mono);

% Scomposizione in sottobande
% x_mono_subband = analysis_fb(x_mono, prototype_dft_filter, I, D);
% Poiché il SFC non è implementaot, il segnale viene replicato sugli L
% canali
% s_subband_base = repmat(x_mono_subband, 1, 1, L);
% Modulazione di fase per decorrelare
s_subband_decorr = phase_modulation_decorrelation(s_subband_base, D, fs);
% Sintesi dei segnati decorrelati in fullband
if (x_test)
    x_speakers = randn(N_recon, L);
else
    x_speakers = zeros(N_recon, L);
    for l = 1:L
        x_speakers(:, l) = synthesis_fb(s_subband_decorr(:, :, l), prototype_dft_filter, I, D);
    end
end
% Generazione dei segnali catturati dai microfoni
y = zeros(N_recon, M);
for m=1:M
    for l=1:L
        % y(:,m) = y(:,m) + filter(H(:, m, l), 1, x);
        y(:,m) = y(:,m) + filter(H(:, m, l), 1, x_speakers(:, l));
    end
end

s_subband = zeros(K_len, I, L);
y_subband = zeros(K_len, I, M);

% Applicazione del banco di analisi a tutti i canali
for l = 1:L
    s_subband(:,:,l) = analysis_fb(x_speakers(:,l), prototype_dft_filter, I, D);
end

for m = 1:M
    y_subband(:,:,m) = analysis_fb(y(:,m), prototype_dft_filter, I, D);
end

% Implementazione dei pesi
sigma_i = zeros(I, 1);
%epsilon_w = 1e-5;

for i = 1:I
    % Calcola la potenza media sull'intero asse temporale (dim 1) e su tutti i canali (dim 3)
    sigma_i(i) = mean(abs(s_subband(:, i, :)).^2, 'all'); 
end

epsilon_w = 0.01 * max(sigma_i);

% Calcolo dei pesi inversamente proporzionali all'energia e normalizzati
w_raw = 1 ./ (sigma_i + epsilon_w);
w_i = w_raw / mean(w_raw); % Normalizzazione a 1

s_i_state = complex(zeros(Ki, I, L));       % Stato per l'impilamento dei campioni
S_i_memory = complex(zeros(L*Ki, P, I));    % Matrice di decorrelazione con P stati passati
S_i_mat = zeros(L*Ki, P);
H_subband = complex(zeros(I, M, L*Ki));     % Matrici dei filtri adattativi: I x M x (L*Ki)
s_ij = zeros(L*Ki, 1);
y_ki = zeros(M, 1);
e_ki = zeros(M, 1);
norm_mis = zeros(1, floor(K_len/100));

u_ij = zeros(Ki,1);

% Calibrazione per il calcolo dei ritardi tramite impulso di test
impulse_len = K + 2*V;
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

for k=1:K_len
    for i=1:I
        s_ij = complex(zeros(L*Ki, 1));
        for l = 1:L
             s_i_state(:, i, l) = [s_subband(k, i, l); s_i_state(1:end-1, i, l)];
             s_ij((l-1)*Ki+1:l*Ki) = s_i_state(:, i, l);
        end

        norm_s = norm(s_ij)^2;
        if norm_s < 1e-6
            % Nei file digitali possono esserci dei silenzi che fanno
            % divergere l'algoritmo. Se la norma è troppo bassa, passare al
            % prossimo campione
            S_i_mat = S_i_memory(:, :, i);
            S_i_memory(:, :, i) = [s_ij, S_i_mat(:, 1:end-1)];
            continue; 
        end

        % In MATLAB invece di inv() è consigliabile usare l'operatore \ che
        % risolve per x il sistema Ax=b

        % A rigore, si deve calcolare a = A^{-1}b con A = (S_j' * S_j) e b
        % = (S_j' * u_ij). Questo sistema può essere portato nella forma
        % Ax=b moltiplicando entrambi i membi per A, e ottenendo quindi
        % (S_j' * S_j) a = (S_j' * u_ij). Per funzionamento dell'operatore
        % \, quindi, a = (S_j' * S_j) \ (S_j'* s_ij).

        % L'aggiunta di reg_ap*eye(P) somma un valore proporzionale all'energia di S_l
        % sulla diagonale di S_j' * S_j, spostando i suoi autovalori dallo
        % zero e evitando instabilità dovuta al fatto che, per segnali correlati, S_j' *
        % S_j è quasi singolare
        if P > 0
            % Estrazione della memoria passata congiunta ((L*Ki) x P)
            S_i_mat = S_i_memory(:, :, i);
            
            % Regolarizzazione della matrice di autocorrelazione
            reg_ap = delta_ap * trace(S_i_mat' * S_i_mat) / P + 1e-6;
            
            % Calcolo dei coefficienti di predizione a_i (P x 1)
            a_i = (S_i_mat' * S_i_mat + reg_ap * eye(P)) \ (S_i_mat' * s_ij);
            
            % Segnale decorrelato congiunto
            u_ij = s_ij - S_i_mat * a_i;
            
            % Aggiornamento della memoria passata
            S_i_memory(:, :, i) = [s_ij, S_i_mat(:, 1:end-1)];
        else
            u_ij = s_ij; % Nessuna decorrelazione
        end
        % Calcolo dell'errore
        y_ki = squeeze(y_subband(k, i, :)); % Estrazione del vettore Mx1
        if M==1, y_ki = y_ki'; end
        H_hat = reshape(H_subband(i,:,:), [M, L*Ki]); % Estrazione della matrice MxL*Ki
        e_ki = y_ki - H_hat * s_ij; % Errore nella sottobanda i-esima
        norm_u = norm(u_ij)^2;
        H_hat = H_hat + mu_h * w_i(i) * (e_ki * u_ij') / (norm_u + delta_h);

        % Aggiornamento della matrice
        H_subband(i, :, :) = H_hat;
    end
    
    % Normalized misalignment calcolato ogni 100 campioni per limitare la
    % complessità
    if (normalized_mis_flag && mod(k, norm_iteration_factor) == 0)
        % Ricostruzione della H fullband:
        % per la ricostruzione corretta si fa passare l'impulso diviso in
        % sottobande attraverso le sottobande della RIR, e poi si
        % ricostruisce la RIR partendo da tali uscite

        % In alternativa, la H fullband può essere ricostruita facendo
        % upsampling del segnale in sottobanda e calcolandone la
        % convoluzione con i filtri di analisi e sintesi. Poiché la
        % convoluzione è lineare, è possibile precalcolare la convoluzione
        % dei filtri di analisi e sintesi e poi calcolare una sola
        % convoluzione per ogni sottobanda


        % Implementazione con filtraggio degli impulsi
        % y_sub_imp = zeros(size(sub_test));
        % for i_sub = 1:I
        %     y_sub_imp(:, i_sub) = filter(H_subband(i_sub, :), 1, sub_test(:, i_sub));
        % end
        % H_recon_raw = synthesis_fb(y_sub_imp, prototype_dft_filter, I, D);

        % Implementazione tramite convoluzione col filtro prototipo
        H_recon_full = complex(zeros(K, M, L));
        for m=1:M
            for l=1:L
                H_sub_ml = zeros(I, Ki);
                H_sub_ml(:, :) = H_subband(:, m, ((l-1)*Ki+1):(l*Ki));
                % H_recon_raw = RIR_reconstruction(prototype_dft_filter, H_sub_ml, I, D);
                H_recon_raw = extract_fullband_RIR(H_sub_ml, prototype_dft_filter, I, D, K);
                H_recon_full(:, m, l) = scale_factor * H_recon_raw(delay_calib+1:delay_calib+K);
            end
        end
        norm_mis(k/norm_iteration_factor) = 20*log10(norm(H(:)-H_recon_full(:), 'fro')/norm(H(:), 'fro'));
    end
end

% Ricostruzione della H fullband
H_recon_full = complex(zeros(K, M, L));
for m=1:M
    for l=1:L
        H_sub_ml = zeros(I, Ki);
        H_sub_ml(:, :) = H_subband(:, m, ((l-1)*Ki+1):(l*Ki));
        % H_recon_raw = RIR_reconstruction(prototype_dft_filter, H_sub_ml, I, D);
        H_recon_raw = extract_fullband_RIR(H_sub_ml, prototype_dft_filter, I, D, K);
        H_recon_full(:, m, l) = scale_factor * H_recon_raw(delay_calib+1:delay_calib+K);
    end
end

% Plot del NM
if (normalized_mis_flag)
    figure;
    plot((1:length(norm_mis))*(norm_iteration_factor/fs), norm_mis);
    xlabel('Iteration Time (s)');
    ylabel('NM (dB)');
    title(sprintf('Normalized Misalignment (MIMO %dx%d)', M, L));
    grid on;
end

% Plot di alcune RIR tra altoparlanti e microfoni
for m_plot=1:M_plotted
    for l_plot=1:L_plotted
        H_sub_ml = zeros(I, Ki);
        H_sub_ml(:, :) = H_subband(:, m_plot, ((l_plot-1)*Ki+1):(l_plot*Ki));
        H_recon_raw = RIR_reconstruction(prototype_dft_filter, H_sub_ml, I, D);
        H_recon_example = scale_factor * H_recon_raw(delay_calib + 1 : delay_calib + K);

        figure;
        plot(H(:, m_plot, l_plot), 'LineWidth', 1.5); 
        hold on;
        plot(H_recon_example, 'LineWidth', 1.5);
        legend(sprintf("Real RIR (Mic %d, Spk %d)", m_plot, l_plot), "Estimated RIR");
        title('RIR Time Domain Comparison'); 
        grid on;
    end
end

% Plot degli spettri delle RIR
figure;
K_fft = length(H_recon_example);
H_true_fft = fft(H(:, m_plot, l_plot));
H_est_fft = fft(H_recon_example);
% Si mantengono solo le frequenze positive
half_len = floor(K_fft/2) + 1;
H_true_half = H_true_fft(1:half_len);
H_est_half = H_est_fft(1:half_len);
% Generazione del vettore delle frequenze fisiche in Hz da 0 a fs/2
freq_vec = (0 : half_len - 1) * (fs / K_fft);
mag_true_dB = 20 * log10(abs(H_true_half));
mag_est_dB = 20 * log10(abs(H_est_half));
plot(freq_vec, mag_true_dB, 'LineWidth', 1.5); 
hold on;
plot(freq_vec, mag_est_dB, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('RIR Frequency Response Comparison');
legend(sprintf("Real RIR (Mic %d, Spk %d)", m_plot, l_plot), "Estimated RIR");
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
    % delay: delay introdotto dai banchi filtri

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
        
        % Modulazione complessa del risultato finale fulllband
        t_i = w_base .* exp(1i * 2 * pi * i * t_total_idx / I);
        
        % Accumulo
        w_full = w_full + real(t_i);
    end
    
    % Compensazione del fattore di scala
    w_full = w_full / D;
end

%% Modulazione di fase per decorrelazione del segnale
function z_subband = phase_modulation_decorrelation(s_subband, D, fs)
    % s_subband: matrice di dimensioni [K_len, I, L]
    % D: fattore di decimazione
    % fs: frequenza di campionamento

    [K_len, I, L] = size(s_subband);
    z_subband = zeros(K_len, I, L);

    % Vettore del tempo decimato w (indice dei campioni in sottobanda)
    w = (0 : K_len - 1)';

    % Costruzione degli angoli alpha come mostrati in Tabella I
    alpha_deg = zeros(I, 1);
    for i = 1:I
        if i <= 4          % u = 0-3
            alpha_deg(i) = 20;
        elseif i == 5      % u = 4
            alpha_deg(i) = 40;
        elseif i == 6      % u = 5
            alpha_deg(i) = 70;
        elseif i == 7      % u = 6
            alpha_deg(i) = 90;
        else               % u >= 7
            alpha_deg(i) = 180;
        end
    end
    % Conversione in radianti
    alpha_rad = alpha_deg * (pi / 180);
    % alpha_rad = pi * ones(I, 1);

    for l = 1:L
        % Frequenza di modulazione f_l
        f_l = 2 * (l - 1);

        for i = 1:I
            alpha_u = alpha_rad(i);

            % Calcolo della fase phi in funzione del tempo w
            phi = alpha_u * sin(2 * pi * f_l * w * D / fs);

            % Modulazione
            z_subband(:, i, l) = s_subband(:, i, l) .* exp(1i * phi);
        end
    end
end

%% Ricostruzione della RIR fullband tramite impulso di test
function H_full = extract_fullband_RIR(H_subband_ml, hp, I, D, K)
    % Ricostruisce la RIR testando il sistema di sottobanda con un impulso.
    % Questo emula la fisica esatta di come i segnali passano nel sistema.
    
    impulse_len = K + 2*length(hp);
    x_test = zeros(impulse_len, 1);
    x_test(1) = 1;
    
    % Passaggio banco filtri analisi
    s_sub = analysis_fb(x_test, hp, I, D);
    
    % Applicazione dei filtri di sottobanda identificati
    y_sub = complex(zeros(size(s_sub)));
    for i = 1:I
        y_sub(:, i) = filter(H_subband_ml(i, :), 1, s_sub(:, i));
    end
    
    % Passaggio banco filtri sintesi
    H_full = synthesis_fb(y_sub, hp, I, D);
end