function plot_pcacmp(EPdata, data, EEG, varargin)
% plot_pcacmp - Generic function to plot PCA components in EEGLAB imported from R or ERP PCA toolkit
% Copyright (c) 2022 Andreas Widmann, University of Leipzig
% Author: Andreas Widmann, widmann@uni-leipzig.de

Arg = struct( varargin{ : } );

% Defaults
if ~isfield( Arg, 'colororder' ) || isempty( Arg.colororder )
    Arg.colororder = get( 0, 'DefaultAxesColorOrder' );
end
if ~isfield( Arg, 'alpha' ) || isempty( Arg.alpha )
    Arg.alpha = 0.4;
end
if ~isfield( Arg, 'factors' ) || isempty( Arg.factors )
    Arg.factors = size( EPdata.data, 5 );
end
if ~isfield( Arg, 'condlabels' ) || isempty( Arg.condlabels )
    Arg.condlabels = cellstr( num2str( (1:size( EPdata.data, 3 ) )' ) );
end
if ~isfield( Arg, 'colormap' ) || isempty( Arg.colormap )
    Arg.colormap = 'jet';
end
nFac = length( Arg.factors );

t = EEG.times / 1000;
nCond = size( EPdata.data, 3 );

% Reconstructed data
for iFac = 1:size( EPdata.data, 5 )
    for iSubj = 1:size( EPdata.data, 4 )
        for iCond = 1:size( EPdata.data, 3 )
            pcaData( :, :, iCond, iSubj, iFac ) = EPdata.data( :, 1, iCond, iSubj, iFac ) * EPdata.facVecT( :, iFac )';
        end
    end
end
gavrPcaData = squeeze( mean( pcaData, 4 ) );

% Peak latency
[ ~, peakLatIdx ] = max( EPdata.facVecT );
peakLat = t( peakLatIdx );

if ~isfield( Arg, 'maplimits' ) || isempty( Arg.maplimits ) || ~isfield( Arg, 'rois' ) || isempty( Arg.rois )
    
    % Peak channel
    [ ~, peakChanIdx ] = max( max( abs( squeeze( mean( EPdata.data, 4 ) ) ), [], 2 ), [], 1 );

    % Plot and topo scaling
    maxabs = ceil( max( max( max( abs( mean( data( peakChanIdx, :, :, : ), 4 ) ), [], 3 ), [], 2 ), [], 1 ) );
    if ~isfield( Arg, 'maplimits' ) || isempty( Arg.maplimits )
        Arg.maplimits = [ -maxabs maxabs ];
    end

end
if ~isfield( Arg, 'levellist' ) || isempty( Arg.levellist )
    Arg.levellist =  linspace( Arg.maplimits( 1 ), Arg.maplimits( 2 ), 21 );
end

figure;
for iFac = 1:length( Arg.factors )
    
    if isfield( Arg, 'rois' ) && ~isempty( Arg.rois )
        chanIdx = [];
        for iChan = 1:length( Arg.rois{ iFac } )
            chanIdx( iChan ) = find( strcmp( Arg.rois{ iFac }{ iChan }, { EEG.chanlocs.labels } ) );
        end
%         { EEG.chanlocs( chanIdx ).labels }
    else
        chanIdx = peakChanIdx( Arg.factors( iFac ) );
    end
        
    for iCond = 1:nCond
        
        hPlot( iFac ) = subplot( nFac, 2, ( iFac - 1 ) * 2 + 1 ); %#ok<AGROW>

        plot_ci( t, squeeze( mean( data( chanIdx, :, iCond, : ), 1 ) )', 'Color', 1 - Arg.alpha + Arg.colororder( iCond, : ) * Arg.alpha, 'LineWidth', 1.2  );
        hold all
        
        tmp = squeeze( mean( gavrPcaData( chanIdx, :, iCond, Arg.factors( iFac ) ), 1 ) )';
        h = plot( t, tmp, 'Color', Arg.colororder( iCond, : ), 'LineWidth', 1.2 );
        set( get( get( h, 'Annotation' ), 'LegendInformation'), 'IconDisplayStyle', 'off' );

        hAx = subplot( nFac, nCond * 2, ( iFac - 1 ) * nCond * 2 + nCond + iCond, 'FontSize', 8 );
        EEG.data = gavrPcaData( :, :, iCond, Arg.factors( iFac ) );
%         pop_plotsserpmap( EEG, 'plot', 'erp', 'type', 'scd', 'lambda', 1e-5, 'items', EEG.times( peakLatIdx( Arg.factors( iFac ) ) ), 'maplimits', Arg.maplimits, 'levelList', Arg.levellist, 'colormap', 'jet', 'axishandle', hAx, 'colorbar', 0 )
%         pop_plotsserpmap( EEG, 'plot', 'erp', 'type', 'sp', 'lambda', 0, 'items', EEG.times( peakLatIdx( Arg.factors( iFac ) ) ), 'maplimits', Arg.maplimits, 'levelList', Arg.levellist, 'colormap', 'jet', 'axishandle', hAx, 'colorbar', 0 )
        topoplot( EEG.data( :, peakLatIdx( Arg.factors( iFac ) ) ), EEG.chanlocs, 'maplimits', Arg.maplimits, 'numcontour', Arg.levellist, 'conv', 'on', 'whitebk', 'on', 'colormap', Arg.colormap );
%         hAx.Title.Visible = 'off';
        hAx.Title.String = Arg.condlabels{ iCond };

    end
    
    title( hPlot( iFac ), sprintf( 'Factor %d: %s, %.3f s', Arg.factors( iFac ), char( join( { EEG.chanlocs( chanIdx ).labels }, ', ' ) ), peakLat( Arg.factors( iFac ) ) ), 'FontSize', 8 )
    xlabel( hPlot( iFac ), 'Time [s]' )
    ylabel( hPlot( iFac ), 'ERP amplitude [µV]' )
        
end

% Show legend
if isfield( Arg, 'legend' ) && strcmp( 'on', Arg.legend)
    legend( hPlot( 1 ), Arg.condlabels, 'Location', 'SouthEast' )
end

if ~isfield( Arg, 'yaxislimits' ) || isempty( 'yaxislimits' )
    Arg.yaxislimits = Arg.maplimits;
end
set( hPlot, 'XLim', [ min( t ) max( t ) ], 'YLim', Arg.yaxislimits, 'XGrid', 'on', 'YGrid', 'on', 'YDir', 'reverse', 'PlotBoxAspectRatio', [ 1.6 1 1 ], 'FontSize', 8 )

end
