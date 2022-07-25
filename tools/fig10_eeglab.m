function fig10_eeglab
%%% STEP 3b: Visual Inspection of Topography %%%
% This script enables visual inspection of the factors and reproduces Fig. 10 in MATLAB/EEGLAB
% Author: Andreas Widmann, widmann@uni-leipzig.de
% Copyright (c) 2022 Andreas Widmann, University of Leipzig
%
% Required steps:
%   (1) Adjust filename for EEGLAB template (or fill required values
%       manually).
%   (2) Adjust filename for individual average file (the same PCA was
%       computed on). Complete dimensions in case necessary.
%   (3) Adjust filename for PCA solution exported from R.
%   (4) Reshape and permute individual average matrix. This step is
%       essential and heavily dependent on the format (row order, cf.
%       Fig. 3 "Observations") of the input data. Target order of matrix
%       dimensions is: chan x time x cond x subj.
%   (5) Reshape and permute PCA scores matrix. This step is essential and
%       heavily dependent on the format (row order, cf. Fig. 3
%       "Observations") of the input data. Target order of matrix
%       dimensions is: chan x time x cond x subj x fac.
%   (6) Select conditions to plot (delete other) and compute difference
%       waves in case necessary.

% Basic figure options
figname = 'adPCA_eeglab-1';

% Optional plot_pcacmp arguments; uncomment and adapt
figOpts = {};
figOpts = [ figOpts, { 'condlabels' }, { { { 'sta', 'nov', 'nov-sta' } } } ]; % Condition labels for legend and topography headings
% figOpts = [ figOpts, { 'condlabels' }, { { { 'nov', 'nov-sta' } } } ]; % Condition labels for legend and topography headings
figOpts = [ figOpts, { 'factors' },  { 1:6 } ]; % Which factors to plot
figOpts = [ figOpts, { 'legend' }, { 'on' } ]; % Show legend (on/off)
% figOpts = [ figOpts, { 'colororder' }, { get( 0, 'DefaultAxesColorOrder' ) } ]; % Custom line colors
% figOpts = [ figOpts, { 'alpha' }, { 0.4 } ]; % Grand-average line transparency; note that MATLAB currently does NOT support alpha adjustment for lines and this feature is implemented by color adjustment
% figOpts = [ figOpts, { 'rois' }, { { { { 'F3', 'Fz', 'F4' }, { 'FC1', 'Cz', 'FC2' }, { 'Fz', 'FC1', 'FC2', 'Cz' }, { 'F3', 'Fz', 'F4' }, { 'Cz' }, { 'F3', 'Fz', 'F4' } } } } ]; % ROIs to plot for each factor; default abs factor peak electrode location
% figOpts = [ figOpts, { 'maplimits' }, { [ -10 10 ] } ]; % Color axis limits for topographies; default abs factor peak amplitude
% figOpts = [ figOpts, { 'yaxislimits' }, { [ -8 8 ] } ]; % Y axis limits for ERP plots; default abs grand-average peak amplitude
try colormap = brewermap(64, 'RdBu'); colormap = flipud( colormap ); catch colormap = 'jet'; end
figOpts = [ figOpts, { 'colormap' }, { colormap } ];

% (1) Load an EEGLAB data file for channel locations and some required parameters
TMP = pop_loadset( 'filepath', '../data', 'filename', 'ad-sta.set' );
EEG = eeg_emptyset;
EEG.srate = TMP.srate;
EEG.xmin = TMP.xmin;
EEG.chanlocs = TMP.chanlocs;
EEG.nbchan = TMP.nbchan;
EEG.pnts = TMP.pnts;
EEG.times = ( ( 0:EEG.pnts - 1 ) / EEG.srate + EEG.xmin ) * 1000;
EEG.xmax = max( EEG.times / 1000 );

% (2) Load individual averages
load( '../data/avrdata.mat', 'data', 'dataIdx' )
data = data( dataIdx( :, 1 ) == 1, : ); % Plot adult group only
dataIdx = dataIdx( dataIdx( :, 1 ) == 1, : ); % Plot adult group only
nCond = 2; % 1 = sta, 2 = dev
nSubj = length( unique( dataIdx( :, 3 ) ) );
nChan = length( unique( dataIdx( :, 4 ) ) );

% (3) Load rotated PCA solution
load( '../results/02bc_rotation_score/rotfit_ad23_geomin0.01.mat' )
nFac = size( loadings, 2 );

% (4) Reshape individual averages
data = reshape( data, nChan, nSubj, nCond, EEG.pnts ); % chan x subj x cond x time; depends on order in exported data matrix!
data = permute( data, [ 1 4 3 2 ] ); % chan x time x cond x subj

% (5) Reshape PCA data
EPdata.data = reshape( scores, nChan, nSubj, nCond, nFac ); % chan x subj x cond x fac; depends on order in exported data matrix!
EPdata.data = permute( EPdata.data, [ 1 5 3 2 4 ] ); % chan x time x cond x subj x fac;
EPdata.facVecT = diag(varSD) * loadings;

% (6) Difference waves
EPdata.data( :, :, 3, :, : ) = EPdata.data( :, :, 2, :, : ) - EPdata.data( :, :, 1, :, : );
data( :, :, 3, : ) = data( :, :, 2, : ) - data( :, :, 1, : );
% EPdata.data( :, :, 1, :, : ) = []; % Delete conditions which should not be plotted

% Plot data
plot_pcacmp( EPdata, data, EEG, figOpts{ : } )

% colorbar
set( gcf, 'PaperPosition', [ 2 1.5 18 26.7 ] )
set( gcf, 'PaperOrientation', 'Portrait' )
print( gcf, '-vector', '-dpdf', sprintf( '../results/03b_topoplot_allFactors/%s.pdf', figname ) )
% close( gcf )

end
