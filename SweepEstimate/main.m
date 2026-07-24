clear all
close all
[orig,Fs] = audioread("sweep_10sec_log.wav");
orig=orig';
rec=audioread("misurazione_1.wav");
mic=rec(:,1)';
% mic=mic./max(abs(mic));
len=Fs*10;
delay=4000;
[h]=IRsweep(delay,len,mic,orig);
plot(h)