clear all; close all;clc

Lobo_spec = data("training1/","Lobo",5);
Ketsa_spec = data("training1/","Ketsa",5);
Lately_spec = data("training1/","Lately",5);

[U,S,V,order,threshold,w,sortband1,sortband2,sortband3] = bands_trainer(Lobo_spec,Ketsa_spec,Lately_spec,40);
plotting(threshold,order,[sortband1;sortband2;sortband3],["Lobo","Ketsa","Lately"]);
%%
Lobo_test = data("testing1/","Lobo",1);
size1 = size(Lobo_test,2);
Ketsa_test = data("testing1/","Ketsa",1);
size2 = size(Ketsa_test,2);
Lately_test = data("testing1/","Lately",1);
size3 = size(Lately_test,2);
hiddenLabel = zeros(1, size1+size2+size3);
hiddenLabel(1:size1) = 1;
hiddenLabel(size1+1:size1+size2) = 2;
hiddenLabel(size1+size2+1:end) = 3;

Test_spec = [Lobo_test Ketsa_test Lately_test];
TestMat = U'*Test_spec;  % PCA projection
pval = w'*TestMat;  % LDA projection
testLabel = zeros(1,length(pval));
for i = 1:length(pval)
    if pval(i) <= threshold(1)
        testLabel(i) = order(1);
    elseif pval(i)<=threshold(2)
        testLabel(i) = order(2);
    else
        testLabel(i) = order(3);
    end
end

diff = ceil(1/3*abs(testLabel - hiddenLabel));
sucRate = 1 - sum(diff)/length(diff);
display(sucRate)
%%
% figure(5)
% subplot(1,3,1)
% histogram(sortband3,30); 
% hold on, plot([threshold(1) threshold(1)],[0 30],'r')
% %set(gca,'Xlim',[-3500 1600],'Ylim',[0 30],'Fontsize',14)
% title('Lately')
% subplot(1,3,2)
% histogram(sortband1,30); 
% hold on, plot([threshold(1) threshold(1)],[0 30],'r')
% plot([threshold(2) threshold(2)],[0 30],'r')
% %set(gca,'Xlim',[-3500 1600],'Ylim',[0 30],'Fontsize',14)
% title('Lobo')
% subplot(1,3,3)
% histogram(sortband2,30); 
% hold on, plot([threshold(2) threshold(2)],[0 30],'r')
% 
% %set(gca,'Xlim',[-3500 1600],'Ylim',[0 30],'Fontsize',14)
% title('Ketsa')

function plotting(threshold,order,sortbands,titles)
   small = sortbands(order(1),:);
   small_title = titles(order(1));
   mid = sortbands(order(2),:);
   mid_title = titles(order(2));
   big = sortbands(order(3),:);
   big_title = titles(order(3));
   
   figure()
   subplot(1,3,1)
   histogram(small,30); 
   hold on
   plot([threshold(1) threshold(1)],[0 30],'r')
   title(small_title)
   subplot(1,3,2)
   histogram(mid,30); 
   hold on
   plot([threshold(1) threshold(1)],[0 30],'r')
   plot([threshold(2) threshold(2)],[0 30],'r')
   title(mid_title)
   subplot(1,3,3)
   histogram(big,30); 
   hold on
   plot([threshold(2) threshold(2)],[0 30],'r')
   title(big_title)
end

function [U,S,V,order,threshold,w,sortband1,sortband2,sortband3] = bands_trainer(b1,b2,b3,feature)
    n1 = size(b1,2); 
    n2 = size(b2,2);
    n3 = size(b3,2);
    n = min([n1,n2,n3]);
    b1 = b1(:,1:n);
    b2 = b2(:,1:n);
    b3 = b3(:,1:n);
    
    [U,S,V] = svd([b1 b2 b3],'econ');
    
    bands = S*V'; % projection onto principal components
    U = U(:,1:feature);
   
    band1 = bands(1:feature,1:n);
    band2 = bands(1:feature,n+1:2*n);
    band3 = bands(1:feature,2*n+1:3*n);
    
    
    m1 = mean(band1,2);
    m2 = mean(band2,2);
    m3 = mean(band3,2);
    m = mean(bands(1:feature,:),2);
    
    Sw = 0; % within class variances
    for k=1:n
        Sw = Sw + (band1(:,k)-m1)*(band1(:,k)-m1)';
    end
    for k=1:n
        Sw = Sw + (band2(:,k)-m2)*(band2(:,k)-m2)';
    end
    for k=1:n
        Sw = Sw + (band3(:,k)-m3)*(band3(:,k)-m3)';
    end    
    
    Sb = n*(m1-m)*(m1-m)'+n*(m2-m)*(m2-m)'+n*(m3-m)*(m3-m)'; % between class
    
    [V2,D] = eig(Sb,Sw); % linear discriminant analysis
    [~,ind] = max(abs(diag(D)));
    w = V2(:,ind); w = w/norm(w,2);
    
    vband1 = w'*band1; 
    vband2 = w'*band2;
    vband3 = w'*band3;
    
    sortband1 = sort(vband1);
    sortband2 = sort(vband2);
    sortband3 = sort(vband3);
    
    sortbands = {sortband1, sortband2, sortband3};
    bands_mean = [mean(sortband1) mean(sortband2) mean(sortband3)];
    bands_mean_sorted = sort(bands_mean);
    order = zeros(1,3);
    for i = 1:3
        order(i) = find(bands_mean == bands_mean_sorted(i));
    end
    
    sortband_small = sortbands{1,order(1)};
    sortband_mid = sortbands{1,order(2)};
    sortband_big = sortbands{1,order(3)};
    
    threshold = zeros(1,2);
   
    t1 = length(sortband_small);
    t2 = 1;
    while sortband_small(t1)>sortband_mid(t2)
        t1 = t1-1;
        t2 = t2+1;
    end
    threshold(1) = (sortband_small(t1)+sortband_mid(t2))/2;
    
    t1 = length(sortband_mid);
    t2 = 1;
    while sortband_mid(t1)>sortband_big(t2)
        t1 = t1-1;
        t2 = t2+1;
    end
    threshold(2) = (sortband_mid(t1)+sortband_big(t2))/2;
end

function spec = data(path, band,number_songs)

    spec = [];   
     for j = 1:number_songs
         [song, Fs] = audioread(path+band+num2str(j)+".mp3");           
         Fs = Fs / 4;
            
         song = song(1:4:end,:);
         monosong = mean(song,2);
         monosong = monosong(find(monosong ,1,'first'):find(monosong,1,'last'));
         %5 second intervals
         pieceLength = Fs*5;
         number = floor(length(monosong) / pieceLength); 
         monosong = monosong(1:(pieceLength * number));
         data = reshape(monosong, [pieceLength, number]);
            
         spec = [spec , data];
     end
     spec = abs(fft(spec));
    
end