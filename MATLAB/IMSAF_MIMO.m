clear;
close all;

L = 2; % numero di altoparlanti (numero di canali)
M = 2; % numero di microfoni
M_plotted = 1;
L_plotted = 2;
I = 16; % numero di sottobande
D = 4; % fattore di decimazione
P = 2; % P = 0: nessuna decorrelazione
% frameSize = 4096;
N = 200000; % lunghezza di x
mu_h = 0.03;
delta_h = 1e-2;  
delta_ap = 1e-1;
K = 2048; % lunghezza delle RIR

normalized_mis_flag = true; % Se impostato a true calcola il NM
norm_iteration_factor = 100;

H = zeros(K, M, L);
file_index = 1;
for m=1:M
    for l=1:L
        fid = fopen("IR\IR1_"+string(file_index)+".f64");
        RIR = fread(fid, 'double');
        fclose(fid);
        H_original_peak = max(abs(RIR(1:K)));
        H(:, m, l) = RIR(1:K) / H_original_peak;
        file_index = file_index + 1;
    end
end

% RIR simulata per test
% H = zeros(8192, 1);
% H(10) = 1; 
% H(20) = -0.5; 
% H(40) = 0.2; 


Ki = ceil(K/D);
H_hat = zeros(M, L*Ki);
x = randn(N, L);
y = zeros(N, M);
for m=1:M
    for l=1:L
        y(:,m) = y(:,m) + filter(H(:, m, l), 1, x(:, l));
    end
end

% V = 128;
% prototype_dft_filter = fir1(V-1, 1/I);
% prototype_dft_filter = fir1(V-1, 1/I, kaiser(V, 12));
% prototype_dft_filter = firceqrip(V-1, 1/I, [0.05 0.03]);
% prototype_dft_filter = prototype_dft_filter / sum(prototype_dft_filter);

% Parametri rcosdesign
beta = 0.5;
span = 128;
prototype_dft_filter = rcosdesign(beta, span, I, "sqrt");
prototype_dft_filter = prototype_dft_filter / sqrt(I);

% Plot della risposta del filtro prototipo
figure;
freqz(prototype_dft_filter, 1, 512);

K_len = floor((N - 1) / D) + 1;
s_subband = zeros(K_len, I, L); 
y_subband = zeros(K_len, I, M);

% Applicazione del banco di analisi a tutti i canali
for l = 1:L
    s_subband(:,:,l) = analysis_fb(x(:,l), prototype_dft_filter, I, D);
end
for m = 1:M
    y_subband(:,:,m) = analysis_fb(y(:,m), prototype_dft_filter, I, D);
end

% Implementazione dei pesi
sigma_i = zeros(I, 1);
epsilon_w = 1e-5;

for i = 1:I
    % Calcola la potenza media sull'intero asse temporale (dim 1) e su tutti i canali (dim 3)
    sigma_i(i) = mean(abs(s_subband(:, i, :)).^2, 'all'); 
end

% Calcolo dei pesi inversamente proporzionali all'energia e normalizzati
w_raw = 1 ./ (sigma_i + epsilon_w);
w_i = w_raw / mean(w_raw); % Normalizzazione a 1

s_i_state = zeros(Ki, I, L);       % Stato per l'impilamento dei campioni
S_i_memory = zeros(L*Ki, P, I);    % Matrice di decorrelazione con P stati passati
H_subband = zeros(I, M, L*Ki);     % Matrici dei filtri adattativi: I x M x (L*Ki)
s_ij = zeros(L*Ki, 1);
y_ki = zeros(M, 1);
e_ki = zeros(M, 1);
norm_mis = zeros(1, floor(K_len/100));

u_ij = zeros(Ki,1);

% Calibrazione per il calcolo dei ritardi tramite impulso di test
if exist("V", "var")
    impulse_len = K + 2*V;
else
    impulse_len = K+span;
end
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
        s_ij = zeros(L*Ki, 1);
        for l = 1:L
             s_i_state(:, i, l) = [s_subband(k, i, l); s_i_state(1:end-1, i, l)];
             s_ij((l-1)*Ki+1:l*Ki) = s_i_state(:, i, l);
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
        if P>0
        u_ij = zeros(L*Ki, 1);
            for l = 1:L
                % Estrazione del vettore di stato temporale del canale
                % l-esimo
                s_l = s_i_state(:, i, l);

                % Aggiornamento della memoria specifica del canale
                S_l_memory = S_i_memory((l-1)*Ki+1:l*Ki, :, i);

                % Coefficienti di predizione
                reg_ap = delta_ap * trace(S_l_memory' * S_l_memory) / P + 1e-6;
                a_l = (S_l_memory' * S_l_memory + reg_ap * eye(P)) \ (S_l_memory' * s_l);

                % Calcolo del segnale decorrelato
                u_ij((l-1)*Ki+1:l*Ki) = s_l - S_l_memory * a_l;

                % Aggiornamento della memoria passata del canale l-esimo
                S_i_memory((l-1)*Ki+1:l*Ki, :, i) = [s_l, S_l_memory(:, 1:end-1)];
            end
        else
            u_ij = s_ij; % Nessuna decorrelazione
        end
        % Calcolo dell'errore
        y_ki = squeeze(y_subband(k, i, :)); % Estrazione del vettore Mx1
        if M==1, y_ki = y_ki'; end
        H_hat = reshape(H_subband(i,:,:), [M, L*Ki]); % Estrazione della matrice MxL*Ki
        e_ki = y_ki - H_hat * s_ij; % Errore nella sottobanda i-esima
        norm_u = norm(u_ij)^2;

        % andrebbero implementati i pesi w_i
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
        H_recon_full = zeros(K, M, L);
        for m=1:M
            for l=1:L
                H_sub_ml = zeros(I, Ki);
                H_sub_ml(:, :) = H_subband(:, m, ((l-1)*Ki+1):(l*Ki));
                H_recon_raw = RIR_reconstruction(prototype_dft_filter, H_sub_ml, I, D);
                H_recon_full(:, m, l) = scale_factor * H_recon_raw(delay_calib+1:delay_calib+K);
            end
        end
        norm_mis(k/100) = 20*log10(norm(H(:)-H_recon_full(:), 'fro')/norm(H(:), 'fro'));
    end
end

% Ricostruzione della H fullband
H_recon_full = zeros(K, M, L);
for m=1:M
    for l=1:L
        H_sub_ml = zeros(I, Ki);
        H_sub_ml(:, :) = H_subband(:, m, ((l-1)*Ki+1):(l*Ki));
        H_recon_raw = RIR_reconstruction(prototype_dft_filter, H_sub_ml, I, D);
        H_recon_full(:, m, l) = scale_factor * H_recon_raw(delay_calib+1:delay_calib+K);
    end
end

% Plot del NM
if (normalized_mis_flag)
    figure;
    plot(norm_mis);
    xlabel(sprintf('Iteration Index (x%d)', norm_iteration_factor));
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