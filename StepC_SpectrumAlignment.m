%Run StepC_SpectrumAlignment.m 
%to align the spectra which are in k-space                   -------------%
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

close all;
clear

fileName = 'Dispcali.mat';
load(fileName);
load('Fringe.mat')
%% reload the calibration parameters obtained from last step
MAmean_first = CalStru.MAmean_first;
MAmean_second = CalStru.MAmean_second;
CArray1 = CalStru.CArray1;
CArray2 = CalStru.CArray2;
x1 = CalStru.x1;
x2 = CalStru.x2;
L1 = CalStru.L1;
L2 = CalStru.L2;
%--------------------------------------------------------------------------
pn = 1024*4; %total pixel number of camera lens
NFFT = 131072;% numbers for FFT
for detpIdx = 1:10

    fringe =  Fringe(:,detpIdx);
    fringe_c1 = (fringe(1:2048));
    fringe_c2 = (fringe(2049:4096));
    
    %% resample into k space
    fringe_c1 = interp1(MAmean_first(L1:L2),fringe_c1(L1:L2),x1);
    fringe_c1(isnan(fringe_c1)) = 0;
    fringe_c2 = interp1(MAmean_second(L1:L2),fringe_c2(L1:L2),x2);
    fringe_c2(isnan(fringe_c2)) = 0;
    
    %% calculate the peak position
    depthProf1 = abs(fft(fringe_c1.*CArray1,NFFT));
    [~,maxd1] = max(depthProf1);
    
    depthProf2 = abs(fft(fringe_c2.*CArray2,NFFT));
    [~,maxd2] = max(depthProf2);
    
    fx = [maxd1-1,maxd1,maxd1+1];
    fy = [depthProf1(maxd1-1),depthProf1(maxd1),depthProf1(maxd1+1)];
    pf = polyfit(fx',fy',2);
    peak1(detpIdx) = -(pf(2)/(2*pf(1)));

    fx = [maxd2-1,maxd2,maxd2+1];
    fy = [depthProf2(maxd2-1),depthProf2(maxd2),depthProf2(maxd2+1)];
    pf = polyfit(fx',fy',2);
    peak2(detpIdx) = -(pf(2)/(2*pf(1)));
%     peak2(detpIdx) = maxd2;
end

match = polyfit(peak2,peak1,1);
ls = (1:length(fringe_c2))';
xf2 = ls./match(1); % find the scaling factor
intxf2 = (min(xf2):1:max(xf2))';
phaseMove = repmat((exp(1i.*2.*pi.*(ls./NFFT).*(match(2)./(match(1))))),1,length(fringe_c2(1,:)));
%% the phaseMove should be very small even negligible, will be negected in 
% next step, you can also count it
%%
dx = 1:NFFT;
% find the itentical wavelength alignment point
for AlineIdx = 1:10
    AlineIdx
    st = 786-3*70-1+141-1-50;  %% you might need to estimate where the iteration should start, or you can iterate from the start to the end in full range
%     for locIdx = st:70:786+3*70
    for locIdx = st:st+200
        fringe =  Fringe(:,AlineIdx);
        fringe_c1 = (fringe(1:2048));
        fringe_c2 = (fringe(2049:4096));
        
        %%resample into k space
        fringe_c1 = interp1(MAmean_first(L1:L2),fringe_c1(L1:L2),x1);
        fringe_c1(isnan(fringe_c1)) = 0;
%         fringe_c1 = fringe_c1.*CArray1;
        fringe_c2 = interp1(MAmean_second(L1:L2),fringe_c2(L1:L2),x2);
        fringe_c2(isnan(fringe_c2)) = 0;
%         fringe_c2 = fringe_c2.*CArray2;
        fringe_c2r = interp1(xf2,fringe_c2,intxf2);% resample
        fringe_c2r(isnan(fringe_c2r)) = 0;
    
        KAline1 = fringe_c1;
        KAline1(length(fringe_c1)+1:length(fringe_c1)+1-locIdx+length(fringe_c2r)-1) = 0;
        KAline2 = zeros(length(fringe_c2r)+100,1);
        KAline2(length(fringe_c1)+1-locIdx:length(fringe_c1)+1-locIdx+length(fringe_c2r)-1) = fringe_c2r;

        KAline = KAline1+KAline2;
        MAG1(locIdx-st+1,AlineIdx) = sum(KAline.^2);

        figure(1); hold on
        plot(KAline1); plot(KAline2);
        close 1;
    end
        fringe1(AlineIdx,:) = fringe_c1;
        fringe2(AlineIdx,:) = fringe_c2;
end
figure(10);
plot(MAG1./repmat(mean(MAG1,1),[201 1])-1)

%%based on this figure, we can find the shift should be the 45, it is also
%%tempatative to chose 44. if we interpolate the fringe into  4 times the
%%length of the fringe, the shift will be more precisely determined.
