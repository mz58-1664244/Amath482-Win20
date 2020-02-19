clear all; close all;clc

%% set up
load handel
v = y';
% plot((1:length(v))/Fs,v);
% xlabel('Time [sec]');
% ylabel('Amplitude');
% title('Signal of Interest, v(n)');

n = length(v);
t = (1:n)/Fs;
v = v(1:n);
L = t(end);
k = 2*pi/L*[0:(n-1)/2 -(n-1)/2:-1];
ks = fftshift(k);
%% Spectrograms for varying Gaussian window width
tslide = 0:0.1:L;
a_vec = [10 1 0.1];

for j = 1:length(a_vec)
    a = a_vec(j);
    vgt_spec = zeros(length(tslide),n);
    for i = 1:length(tslide)
        g = exp(-a*(t-tslide(i)).^2);
        vg = g.*v;
        vgt = fft(vg);
        vgt_spec(i,:) = fftshift(abs(vgt));
    end
    subplot(1,3,j)    
    pcolor(tslide,ks,vgt_spec.'),     
    shading interp     
    title(['a = ',num2str(a)])   
    xlabel('Time'), ylabel('Frequency')
        
    colormap(hot) 
end
%% Spectrograms for varying Gaussian window translation

steps = [0.01 0.1 1];
for j = 1:length(steps)
    step = steps(j);
    tslide = 0:step:9;
    vgt_spec = zeros(length(tslide),n);    
    for i = 1:length(tslide)
        g = exp(-(t-tslide(i)).^2);
        vg = g.*v;
        vgt = fft(vg);
        vgt_spec(i,:) = fftshift(abs(vgt));
    end
    subplot(1,3,j)    
    pcolor(tslide,ks,vgt_spec.'),     
    shading interp     
    title(['step = ',num2str(step)])   
    xlabel('Time'), ylabel('Frequency')
        
    colormap(hot) 
end
%% Spectrograms for different Gabor windows
tslide = 0:0.1:L;
name = ["Gaussian","Mexican","Shannon"];
for j = 1:3
    vgt_spec = zeros(length(tslide),n); 
    for i = 1:length(tslide)
        if (j == 1)
            g = exp(-(t-tslide(i)).^2);        
        elseif (j == 2)
            g = (1-(t-tslide(i)).^2).*exp(-(t-tslide(i)).^2);
        else
            g = abs(t-tslide(i)) <= 0.5;
        end
        vg = g.*v;
        vgt = fft(vg);
        vgt_spec(i,:) = fftshift(abs(vgt));
    end
    subplot(1,3,j)    
    pcolor(tslide,ks,vgt_spec.'),     
    shading interp     
    title(name(j))   
    xlabel('Time'), ylabel('Frequency')
        
    colormap(hot)
end
        


