clear all
close all
clc
%% DATA INPUT: MASW wavefield 
Filename = 'BH01_10m_offset_F.dat'; % Site-1,Forward Shot
HeaderLines = 7;
fs = 8000;  %Sampling frequency (Hz)
N = 24;    % Number of Geophones
x1 = 10; % Offset Distance (m)
dx = 1; % Receiver Spacing
Direction = 'forward'; %Geophone taken in forward direction from 1 to 24 

[u,t,Tmax,L,x] = MASWaves_read_data(Filename,HeaderLines,fs,N,dx,x1,Direction);

%% PLOT DATA 
du = 1/20; % wavefield scale
FigWidth = 10; % cm
FigHeight = 8; % cm
FigFontSize = 14; % pt
figure
MASWaves_plot_data(u,N,dx,x1,L,t,Tmax,du,FigWidth,FigHeight,FigFontSize)

%% DISPERSION IMAGE GENERATION
cT_min = 100; % m/s
cT_max = 500; % m/s
delta_cT = 0.5; % m/s

[f,c,A] = MASWaves_dispersion_imaging(u,N,x,fs,cT_min,cT_max,delta_cT);

%% VELOCITY SPECTRA 2D
resolution = 100; 
fmin = 0; % Hz
fmax = 50; % Hz
FigWidth = 10; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
[fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_2D(f,c,A,fmin,fmax,...
    resolution,FigWidth,FigHeight,FigFontSize);

%% VELOCITY SPECTRA 3D
fmin = 0; % Hz
fmax = 50; % Hz
FigWidth = 10; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
[fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_3D(f,c,A,fmin,fmax,...
    FigWidth,FigHeight,FigFontSize);

%% SELECT POINTS 
f_receivers = 4.5; % Hz
select = 'numbers'; % After run enter "1:14" to select the point in fundamental mode
up_low_boundary = 'no'; 
p = 95; % Percentage
[f_curve0,c_curve0,lambda_curve0,...
    f_curve0_up,c_curve0_up,lambda_curve0_up,...
    f_curve0_low,c_curve0_low,lambda_curve0_low] = ...
    MASWaves_extract_dispersion_curve(f,c,A,fmin,fmax,f_receivers,...
    select,up_low_boundary,p);

%% DISPERSION CURVE GENERATION
FigWidth = 9; % cm
FigHeight = 6; % cm
FigFontSize = 8; % pt
type = 'f_c';
up_low_boundary = 'yes';
figure
MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
     f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,c_curve0_low,...
     lambda_curve0_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize)

FigWidth = 7; % cm
FigHeight = 9; % cm
FigFontSize = 12; % pt
type = 'c_lambda';
up_low_boundary = 'yes';
figure
MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
     f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,c_curve0_low,...
     lambda_curve0_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize)
 

%% Inversion Analysis to estimate shear wave velocity profile
c_test_min = 0; % m/s
c_test_max = 700; % m/s
delta_c_test = 0.5; % m/s
c_test = c_test_min:delta_c_test:c_test_max; % m/s

% Layer parameters
n = 6; 
beta =[152	152	227.8	310	181.5	310 310]; % Optimized Vs values using TLBO algorithm 
h =[0.3	1.8	1.0	1.8	3.0	3.6 1]; % Optimized thickness using TLBO algorithm
alpha =[1000 1000 1000 1000 1000 1000 1000]; % Assumed compressional wave velocity (m/s)
rho = [2000 2000 2000 2000 2000 2000 2000]; % Assumed soil density (kg/m^3)

up_low_boundary = 'yes';
[c_t,lambda_t] = MASWaves_theoretical_dispersion_curve...
    (c_test,lambda_curve0,h,alpha,beta,rho,n);

up_low_boundary = 'yes';
FigWidth = 8; % cm
FigHeight = 10; % cm
FigFontSize = 12; % pt
f_curvet = f_curve0';
figure
MASWaves_plot_theor_exp_dispersion_curves(c_t,lambda_t,...
    c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
    c_curve0_low,lambda_curve0_low,up_low_boundary,...
    FigWidth,FigHeight,FigFontSize)

e = MASWaves_misfit(c_t,c_curve0);

up_low_boundary = 'yes';
FigWidth = 16; % cm
FigHeight = 10; % cm
FigFontSize = 12; % pt
figure
MASWaves_plot_inversion_results_one_iteation(c_t,f_curvet,...
    c_curve0,f_curve0,c_curve0_up,f_curve0_up,c_curve0_low,...
    f_curve0_low,n,beta,h,e,up_low_boundary,FigWidth,FigHeight,FigFontSize)