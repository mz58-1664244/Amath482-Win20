clear all; close all; clc; 
%% Case1
load('cam1_1.mat') 
load('cam2_1.mat') 
load('cam3_1.mat')
data1_1= dataExtraction(vidFrames1_1, 170:430, 300:400, 250);
data2_1= dataExtraction(vidFrames2_1, 100:400, 225:355, 250);
data3_1= dataExtraction(vidFrames3_1, 200:350, 235:480, 247);
syncData1 = syncData(data1_1, data2_1, data3_1);
PCAplot(syncData1);

%% Case2
load('cam1_2.mat') 
load('cam2_2.mat') 
load('cam3_2.mat')
data1_2= dataExtraction(vidFrames1_2, 170:430, 300:450, 250);
data2_2= dataExtraction(vidFrames2_2, 50:475, 165:450, 249);
data3_2= dataExtraction(vidFrames3_2, 180:350, 235:500, 246);
syncData2 = syncData(data1_2, data2_2, data3_2);
PCAplot(syncData2);

%% Case3
load('cam1_3.mat') 
load('cam2_3.mat') 
load('cam3_3.mat')
data1_3= dataExtraction(vidFrames1_3, 225:450, 275:450, 250);
data2_3= dataExtraction(vidFrames2_3, 100:425, 165:415, 249);
data3_3= dataExtraction(vidFrames3_3, 160:365, 235:495, 246);
syncData3 = syncData(data1_3, data2_3, data3_3);
PCAplot(syncData3);

%% Case4
load('cam1_4.mat') 
load('cam2_4.mat') 
load('cam3_4.mat')
data1_4= dataExtraction(vidFrames1_4, 225:450, 275:470, 245);
data2_4= dataExtraction(vidFrames2_4, 50:400, 165:425, 249);
data3_4= dataExtraction(vidFrames3_4, 100:300, 270:505, 230);
syncData4 = syncData(data1_4, data2_4, data3_4);
PCAplot(syncData4);
%% functions
function data = dataExtraction(vidFrames, filter_y, filter_x, thresh)
    %Play video
    %implay(vidFrames);
    
    filter = zeros(480,640);
    filter(filter_y, filter_x) = 1;
    data = [];
    numFrames = size(vidFrames,4);
    %figure()
    for j=1:numFrames
        X=vidFrames(:,:,:,j);
        Xgrey = rgb2gray(X);
        Xgrey = double(Xgrey);
        Xf = Xgrey.*filter;
        threshold = Xf > thresh;
        indeces = find(threshold);
        [y, x] = ind2sub(size(threshold),indeces);
        data = [data, [mean(x); mean(y)]];        
%       subplot(1,2,1)
%       imshow(uint8((threshold * 255))); 
%       drawnow 
%       subplot(1,2,2)
%       imshow(uint8(Xf)); drawnow
    end
end

function syncData = syncData(data1, data2, data3)
    [minimum, startIndex]= min(data1(2,1:20));
    data1 = data1(:, startIndex:end);
    [minimum, startIndex]= min(data2(2,1:20));
    data2 = data2(:, startIndex:end);
    [minimum, startIndex]= min(data3(1,1:20));
    data3 = data3(:, startIndex:end);
    len = min([length(data1) length(data2) length(data3)]);
    data1 = data1(:, 1:len);
    data2 = data2(:, 1:len);
    data3 = data3(:, 1:len);
    syncData = [data1;data2;data3];
end

function PCAplot(syncData)
    [m,n]=size(syncData);  
    mn=mean(syncData,2); 
    syncData =syncData-repmat(mn,1,n); 
    [u,s,v]=svd(syncData/sqrt(n-1), 'econ'); 
    lambda=diag(s).^2;
    project= u'*syncData ; 
    sig=diag(s);

    figure(1)
    clf;
   
    % Plot energy
    subplot(1,2,1)
    semilogy(lambda/sum(lambda), 'ro', 'linewidth', 2)
    xlabel('PC number')
    ylabel('energy percent(log scale)')

    % Plot Displacement in principal component direction
    subplot(1,2,2)
    hold on
    plot(1:length(project), project(1:3,:),'linewidth', 1.5)
    xlabel('Time')
    ylabel('Displacement')
    xlim([0, length(project)])
    legend('PC 1','PC 2', 'PC 3')
end