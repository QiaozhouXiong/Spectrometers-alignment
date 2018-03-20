%Run StepB_RSK_DispC.m 
%to routine calibration to transform spectrum into k-space and compensate dispersion-%
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
%-------------------------------------------------------------------------%
%% this part is software method to transform the spectrum from wavelength 
%space into k-space, and doin the dispersion calibration, for more detais
%in this part, pls refer to X. Yu, X. Liu, S. Chen, Y. Luo, X. Wang, 
%and L. Liu, "High-resolution extended source optical coherence tomography," Opt. Express 23, 26399 (2015).
%%------------------------------------------------------------------------%%



%%
close all;
clear;
tic
%% Seting the parameters for dispersion calibration
NFFT = 4096*32;
t = 1:2048;
T1 = [600 600];
T2 = [1200 1400];
startPo = [262 50];
endPo = [1885 2000];
L1 = startPo;
L2 = endPo;
%%-------------------------------------------------------------------------
pn = 1024*4; %total pixel number of camera lens

fileName = 'Dispcali.mat';
% load the background of the fringes
s = (double(importdata('bgn.mat')));

%% read the interference
fileMark = dir('*.mat');
Frame_No = 10; 
Fringe = zeros(pn,Frame_No);
for index = 1:Frame_No
    filemat = sprintf('f%d.mat',index);
    Temp = importdata(filemat);
    Fringe(:,index) = Temp(:,1);
end

ss = repmat(s,[1 Frame_No]);
Fringe = Fringe - ss;

Fringe(1:startPo(1),:) = 0;
Fringe(endPo(1):2048+startPo(2),:) = 0;
Fringe(2048+endPo(2):end,:) = 0;
save('Fringe','Fringe')
MA_first = zeros(2048,Frame_No/2);
MA_second = zeros(2048,Frame_No/2);

%% transform the spectrum into evenly sampled in k-space by hilbert transform.
for index = 1:Frame_No/2

    %%------------------pure fringe--------------------------------------
    First_fringe =  Fringe(1:2048,index);
    Second_fringe =  Fringe(2049:end,index);
    
    y_first1 = First_fringe;
    Hy_first1 = hilbert(y_first1);
    HyAng_first1 = unwrap(angle(Hy_first1));
    HyAng_frist_net1 = HyAng_first1 - HyAng_first1(T1(1));
    
    y_second1 = Second_fringe;
    Hy_second1 = hilbert(y_second1);
    HyAng_second1 = unwrap(angle(Hy_second1));
    HyAng_second_net1 = HyAng_second1 - HyAng_second1(T1(2));
    
    
    First_fringe =  Fringe(1:2048,index+Frame_No/2);
    Second_fringe =  Fringe(2049:end,index+Frame_No/2);
    
    y_first2 = First_fringe;
    Hy_first2 = hilbert(y_first2);
    HyAng_first2 = unwrap(angle(Hy_first2));
    HyAng_frist_net2 = HyAng_first2 - HyAng_first2(T1(1));
    
    y_second2 = Second_fringe;
    Hy_second2 = hilbert(y_second2);
    HyAng_second2 = unwrap(angle(Hy_second2));
    HyAng_second_net2 = HyAng_second2 - HyAng_second2(T1(2));
    
    MappingAng_frist = HyAng_frist_net2 - HyAng_frist_net1;
    MappingAng_frist = MappingAng_frist/MappingAng_frist(T2(1));
    MappingAng_second = HyAng_second_net2 - HyAng_second_net1;
    MappingAng_second = MappingAng_second/MappingAng_second(T2(2));
    
    MA_first(:,index) = MappingAng_frist;
    MA_second(:,index) = MappingAng_second;
    
end
clear AngMap_frist;
clear AngMap_second;
%%
MAmean_first = mean(MA_first,2);
MAmean_second = mean(MA_second,2);

x1 = (linspace(MAmean_first(L1(1)),MAmean_first(L2(1)),L2(1)-L1(1)+1))';
x2 = (linspace(MAmean_second(L1(2)),MAmean_second(L2(2)),L2(2)-L1(2)+1))';


%% Calibrate the dispersion
CPhaseArray1 = zeros(L2(1)-L1(1)+1,Frame_No);
CPhaseArray2 = zeros(L2(2)-L1(2)+1,Frame_No);

CP1 = zeros(L2(1)-L1(1)+1,Frame_No);
CP2 = zeros(L2(2)-L1(2)+1,Frame_No);
for index = 1:Frame_No
    
    First_fringe =  Fringe(1:2048,index);
    Second_fringe =  Fringe(2049:end,index);
    
    y_first = First_fringe;
    yC_first = interp1(MAmean_first(L1(1):L2(1)),y_first(L1(1):L2(1)),x1);
    
    y_second = Second_fringe;
    yC_second = interp1(MAmean_second(L1(2):L2(2)),y_second(L1(2):L2(2)),x2);
    
    yCAng_first = unwrap(angle(hilbert(yC_first)));
    p1 = polyfit((0:(length(yCAng_first)-1))',yCAng_first,1);
    CPhase1 = yCAng_first - polyval(p1,(0:(length(yCAng_first)-1))');
    
    yCAng_second = unwrap(angle(hilbert(yC_second)));
    p2 = polyfit((0:(length(yCAng_second)-1))',yCAng_second,1);
    CPhase2 = yCAng_second - polyval(p2,(0:(length(yCAng_second)-1))');
    
    CP1(:,index) = CPhase1;
    CP2(:,index) = CPhase2;
end
%%

%% set the limit of dispersion calibration
CPhase1 = mean(CP1,2);
CPhase1(CPhase1>5) = 5;
CPhase1(CPhase1<-5) = -5;
CArray1 = exp(-(1i.*CPhase1));

CPhase2 = mean(CP2(:,1:6),2);
CPhase2(CPhase2>5) = 5;
CPhase2(CPhase2<-5) = -5;
CArray2 = exp(-(1i.*CPhase2));

%% save the data into one struct format
CalStru.MAmean_first = MAmean_first;
CalStru.MAmean_second = MAmean_second;
CalStru.x1 = x1;
CalStru.x2 = x2;
CalStru.CArray1 = CArray1;
CalStru.CArray2 = CArray2;
CalStru.L2 = L2;
CalStru.L1 = L1;
CalStru.NFFT = 2048;
save(fileName,'CalStru');
toc;

