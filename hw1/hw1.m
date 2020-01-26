clear; close all; clc;
%% load data and set up the domain
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);
% for j=1:20
%    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
%    close all, isosurface(X,Y,Z,abs(Un),0.4)
%    axis([-20 20 -20 20 -20 20]), grid on, drawnow
%    pause(1)
% end
%% Average spectrums and Find Center Frequency
avgUnt = zeros(64,64,64);
for j=1:20
   Un(:,:,:)=reshape(Undata(j,:),n,n,n);
   Unt = fftn(Un);
   avgUnt(:,:,:) = avgUnt + Unt;
end
avgUnt = abs(fftshift(avgUnt))/20;
[maximum,index] = max(avgUnt(:));
xc = Kx(index)
yc = Ky(index)
zc = Kz(index)

%% Filter Data to get the locations of marble
places = zeros(20,3);

filter = exp(-1.0*fftshift(Kx-xc).^2).* ...
    exp(-1.0*fftshift(Ky-yc).^2) .* ...
    exp(-1.0*fftshift(Kz-zc).^2);

for j = 1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n); 
    Unt = fftn(Un);
    unft = filter.* Unt;    
    unf = ifftn(unft);
    unf = abs(unf);
    [maximum,index] = max(unf(:));
    places(j,:) = [X(index) Y(index) Z(index)];
    %isosurface(X,Y,Z,unf/maximum,0.75)
    %axis([-20 20 -20 20 -20 20]), grid on, drawnow
end
figure()
plot3(places(:,1), places(:,2), places(:,3), 'ro-')
axis([-L L -L L -L L]), xlabel('x'); ylabel('y');zlabel('z');
grid on, drawnow

fplace = places(20,:)