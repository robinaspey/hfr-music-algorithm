%
%   tsreadv13.m
%   Plots filtered data 
%   Needs work to fix FFT of filtered data on bottom row.
%
% Clear everything to start new run
clear all; close all; clc; delete(gcf); clf;

% Initialise variables for filter.
fs          = 62e6;
order       = 10;
fcutlow     = 15.5e6;
fcuthigh    = 16.5e6;

[xb,xa]     = butter(order,[fcutlow,fcuthigh]/(fs/2), 'bandpass');
% close all;
                                  
dirinfo = @getFiles;
% disp('Reading time series file sequence: ');
% dinfo = dir('d:\!data1\*.ts');
% for K = 1 : length(dinfo)
%        thisfilename = dinfo(K).name;
%        % thisdata = load(thisfilename);
%        fprintf('File: #%d, "%s" file found\n', K, thisfilename);
% end;
% pause(1);

formatString = '%f<gtag>%f<atag>%f<indx>%f<scal>%f<alvl>%f' ;
stringy      = 'gtag'                                       ;
fid          =  fopen('d:\\!data0\\test-ts.ts')             ;
n            =  0                                           ;  

if  ( fid <= 0 ) 
    disp('Could not open file; exiting nothing to do !');
    return;
end;

tic;

[A,n]    = fread(fid, 'uchar');
C        = transpose(A)       ;

% Count the number of samples based on the tag locations

m1       = strfind(C, 'gtag') ;
m2       = strfind(C, 'atag') ;
m3       = strfind(C, 'indx') ;
m4       = strfind(C, 'scal') ;
m5       = strfind(C, 'alvl') ;

% Now locate the sizes of the indices arrays of the Key tags

mm1      = max(size(m1));
mm2      = max(size(m2));
mm3      = max(size(m3));
mm4      = max(size(m4));
mm5      = max(size(m5));       
sz       = max(size(m1));
sarr     = max(size(A));

str      = sprintf('Matches:%u, Blocks:%u', sz, n );

disp(str);
toc;

fclose(fid);
% clear all;
% now read in 256 blocks of 24660 bytes or 6312960 bytes total.
% We now have the arrays m1,m2,m3,m4 which have the indices of 
% the tags inside the binary data array C.
disp("-- Debug #1 --");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default maximum is 10
% THIS LOOP CHECK THE DATA FILE BY PRINTING OUT THE VALUES OF THE
% INDEX KEYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ndx=1:5
    % Print indices for key tags inside each Time Series (.ts) 
    % file into arrays to assist in locating the start of data using 
    % the alvl key which is just before the start of data blocks which
    % are 8192 samples long... This seems a little short.
    k1   = m1(ndx);
    k2   = m2(ndx);
    k3   = m3(ndx);
    k4   = m4(ndx);
    k5   = m5(ndx);
    % str1=sprintf('Key=%03d:gtag=%06u,atag=%06u,indx=%06u,..
    %     scal=%06u,alvl=%06u;%06ld,%06ld,%06ld,%06ld,%06ld',...
    %               ndx,k1,k2,k3,k4,k5,k1v,k2v,k3v,k4v,k5v);
    % Now print some values out for checking that the files are being 
    % read consistently in terms of the index values.
    str1  = sprintf('Key=%03d:gtag=%06u,atag=%06u,indx=%06u,scal=%06u,alvl=%06u', ndx,k1,k2,k3,k4,k5);
    disp( str1 );
    pause(0.1);
end
% Now extract data blocks starting from offset of first alvl at 397
% in sets of 3 with length of 8192 uint8 or 4096 of uint16.
% First initialise data arrays for L1, L2, L3 loops 
% Find the size of the array ......................................
maxval   = max(size(m5));
% define holders for data for loop1, loop2 and whip
loop1    = zeros(8192);       
loop2    = zeros(8192);
loop3    = zeros(8192);

ndxL1    = 0;
ndxL2    = 0;
ndxL3    = 0;
%
%   Set the sampling frequency
%
Fs       = 62e6;                  % default 60MHz
T        = 1/Fs;                  % Sampling period
L        = 4096;                  % 256-8192 length of sampled data
t        = (0:L-1)*T;             %
tic;                              % start clock for time to read all data                     
fig1     =  figure('name','Simulation Plot Window');
set(gcf, 'Position', [100, 100, 1200, 1000])

%   CHANGE: defined xdata16 now as int16
%
xdata16  = int16(4096);
maxslice = 8192/256;

while (ndx<300)
    if (ndx == 0) 
        ndx = 1;
    end;
    
    disp(ndx);
    pause(1);
    
    ds       = m5(ndx+0);   loop1    = C(ds:ds+8191);
    ds       = m5(ndx+1);   loop2    = C(ds:ds+8191);
    ds       = m5(ndx+2);   loop3    = C(ds:ds+8191); % whip antenna
    
    ndx = ndx + 3;
    
    j=1, i=1;

    subplot(4, 3, 1);
    plot(loop1);
    title('P1.DF Loop 1');
    xlabel('sample');
    subplot(4, 3, 2);
    plot(loop2);
    title('P2.DF Loop 2');
    xlabel('sample');
    subplot(4, 3, 3);
    plot(loop3);
    title('P3.Whip');
    xlabel('sample');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now we are calculating the FFT for loop1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L  = 256;
    k  = 0;
    for (k=1:31)
        % Move through buffers for 3 channels of the time series
        % data, the aim is to look for decrease in FFT signal amplitude
        % representative of sea echo for near 16MHz Bragg.
        str=sprintf('Value of k is: %u', k);
        disp(str);
        
        pause(0.1);
        offs          = k*256;
        X1            = loop1(offs:offs+L);
        Y1            = fft(X1,L);
        P1            = abs(Y1/L);
        fft1          = P1(1:L/2+1);
    
        fft1(2:end-1) = 2*fft1(2:end-1);
        scalex        = Fs*(0:(L/2))/L;
        len           = max(size(scalex));
    
        subplot(4, 3, 4);
        plot(scalex(4:len)/1e6, fft1(4:len));
        ylim([0 100]);
        str=sprintf('P4.Loop1:%u (offset)', offs);
        title(str);
        xlabel('Freq/MHz');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Now we are calculating the FFT for loop2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X1            = loop2(offs:offs+L);
        Y1            = fft(X1,L);
        P1            = abs(Y1/L);
        fft1          = P1(1:L/2+1);
    
        fft1(2:end-1) = 2*fft1(2:end-1);
        scalex        = Fs*(0:(L/2))/L;
        len           = max(size(scalex));
    
        subplot(4, 3, 5);

        plot(scalex(4:len)/1e6, fft1(4:len));
        ylim([0 100]);
        str=sprintf('P5.Loop2:%u (offset)', offs);
        title(str);
        xlabel('Freq/MHz');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Now we are calculating the FFT for loop3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X1            = loop3(offs:offs+L);
        Y1            = fft(X1, L);
        P1            = abs(Y1/L);
        fft1          = P1(1:L/2+1);
    
        fft1(2:end-1) = 2*fft1(2:end-1);
        scalex        = Fs*(0:(L/2))/L;
        len           = max(size(scalex));
    
        subplot(4, 3, 6);
        plot(scalex(4:len)/1e6, fft1(4:len));
        ylim([0 100]);
        str = sprintf('P6.Whip:%u (offset)', offs);
        title(str);
        
        xlabel('Freq/MHz');
        
        x1        = filter(xb,xa,loop1);
        x2        = filter(xb,xa,loop2);
        x3        = filter(xb,xa,loop3);
        
        subplot(4, 3, 7);
        plot(x1);
        xlabel('sample');
        title('P7.BP Filtered')
        ylim([-100 100]);
        subplot(4, 3, 8);
        plot(x2);
        xlabel('sample');
        title('P8.BP Filtered');
        ylim([-100 100]);
        subplot(4, 3, 9);
        plot(x3);
        xlabel('sample');
        title('P9.BP Filtered');
        ylim([-100 100]);
        
        % L=1024;
        
        subplot(4, 3, 10);
        
        n = 8192;
        
        Y = fft(x1,n);
        % Define the frequency domain and plot the unique frequencies.

        f = 0.1*Fs*(0:(n/2))/n;
        P = abs(Y/n);

        plot(f,P(1:n/2+1)) 
        xlim([1.5e6 1.7e6]);
        ylim([0 12]);
        title('P10.FFT 1')
        xlabel('Freq/MHz')
        ylabel('P(f)')
        
        subplot(4, 3, 11);
        
        Y = fft(x2,n);
        % Define the frequency domain and plot the unique frequencies.

        f = 0.1*Fs*(0:(n/2))/n;
        P = abs(Y/n);

        plot(f,P(1:n/2+1)) 
        xlim([1.5e6 1.7e6]);
        ylim([0 12]);
        title('P11.FFT 2')
        xlabel('Freq/MHz')
        ylabel('P(f)')
        
        subplot(4, 3, 12);
        
        Y = fft(x3,n);
        % Define the frequency domain and plot the unique frequencies.

        f = 0.1*Fs*(0:(n/2))/n;
        P = abs(Y/n);
        xlim([1e6 5e6]);
        plot(f,P(1:n/2+1)) 
        xlim([1.5e6 1.7e6]);
        ylim([0 12]);
        title('P12.FFT 3')
        xlabel('Freq/MHz')
        ylabel('P(f)')

    end
    % figure('name','Hermitian Matrix');
    % hermityy = (1/(Y*htranspose(Y)));
    % str=sprintf(' Hermitian size: ');
    % disp(str);
    % size(hermityy);
end 
toc;

function dinfo = getfiles()
    disp('Reading time series file sequence: ');
    dinfo = dir('d:\!data1\*.ts');

    for K = 1 : length(dinfo)
        thisfilename = dinfo(K).name;
        % thisdata = load(thisfilename);
        fprintf('File: #%d, "%s" file found\n', K, thisfilename);
    end
    pause(1);
end




