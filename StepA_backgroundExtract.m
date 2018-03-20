%Run StepA_backgroundExtract.m 
%to extract background offset within the interference signals-------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

%Author: Qiaozhou Xiong(E150022@e.ntu.edu.sg)
%Affiliation: School of EEE, Nanyang Technological Univerisity
%Lastest revision: Feb 16 2018 / Last Comment revision Feb 16 2018

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%for more details, pls refer to our publication, citation appreciated but not required.
% Q.Xiong,et al, A generic method to co-register two spectrometers in------ 
% spectral domain optical coherence tomography 
%%

clc
clear
%% rread

pn = 1024*4; %total pixel number of camera lens
c = 512*2;     %number of line scans in a frame

fileMark = dir('*.mat');

bgnsum = zeros(pn,1);
frameNo = length(fileMark); 
for index = 1:frameNo
    B1 = importdata(fileMark(index).name);
    bgn = mean(B1,2);
    figure(1);
    plot(bgn);
    drawnow;
    bgnsum = bgn+bgnsum;
end
bgn = bgnsum/(frameNo-1+1);
%%

save('bgn','bgn');
