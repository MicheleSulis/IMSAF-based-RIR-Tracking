function [h]=IRsweep(delay,len,mic,orig)

%%% Funzione che calcola la risposta impulsiva con la tecnica dello sweep
%%% delay è il ritardo del segnale registrato
%%% len è la lunghezza del segnale
%%% mic è il segnale microfonico
%%% orig è il segnale inviato (sweep)

orig=orig(1:len);
    
data_FFT=fft(orig);
data_rec_FFT=fft(mic(1,1+delay:len+delay));
h=ifft(data_rec_FFT./data_FFT);
