clear all; close all;clc
[y,Fs] = audioread('music1.wav');
tr_piano=length(y)/Fs;  % record time in seconds

n = length(y);
t = (1:n)/Fs;
y = y(1:n)';
L = t(end);
k = 1/L*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

a = 25;
tslide=0:0.1:L;
ygt_spec = zeros(length(tslide),n);
for j=1:length(tslide)    
    g=exp(-a*(t-tslide(j)).^2);     
    yg=g.*y;     
    ygt=fft(yg);     
    ygt_spec(j,:) = fftshift(abs(ygt)); % We don't want to scale it
end

figure(1)
pcolor(tslide,ks,ygt_spec.'),
shading interp 
xlabel('Time（s)'), ylabel('Frequency(Hz)')
set(gca,'Ylim',[-700 700],'Fontsize',16) 
colormap(hot)
%plot((1:length(y))/Fs,y);
%xlabel('Time [sec]'); ylabel('Amplitude');
%title('Mary had a little lamb (piano)');
%p8 = audioplayer(y,Fs); playblocking(p8);

%%
%figure(2)
[y,Fs] = audioread('music2.wav');
tr_rec=length(y)/Fs;  % record time in seconds
n = length(y);
t = (1:n)/Fs;
y = y(1:n)';
L = t(end);
k = 1/L*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

a = 25;
tslide=0:0.1:L;
ygt_spec = zeros(length(tslide),n);
for j=1:length(tslide)    
    g=exp(-a*(t-tslide(j)).^2);     
    yg=g.*y;     
    ygt=fft(yg);     
    ygt_spec(j,:) = fftshift(abs(ygt)); % We don't want to scale it
end

figure(1)
pcolor(tslide,ks,ygt_spec.'),
shading interp 
xlabel('Time（s)'), ylabel('Frequency(Hz)')
set(gca,'Ylim',[-2000 2000],'Fontsize',16) 
colormap(hot)
%plot((1:length(y))/Fs,y);
%xlabel('Time [sec]'); ylabel('Amplitude');
%title('Mary had a little lamb (recorder)');
%p8 = audioplayer(y,Fs); playblocking(p8);