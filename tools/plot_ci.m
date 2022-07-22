function h = plot_ci( x, y, varargin )
% plot_ci - Plot with confidence interval
% Copyright (c) 2022 Andreas Widmann, University of Leipzig
% Author: Andreas Widmann, widmann@uni-leipzig.de

fieldArray =  lower( varargin( 1:2:end - 1 ) );
argArray =  varargin( 2:2:end );

Arg = cell2struct( argArray, fieldArray, 2 );

if ~isfield( Arg, 'dim' ) || isempty( Arg.dim )
    Arg.dim = 1;
end
if ~isfield( Arg, 'facealpha' ) || isempty( Arg.facealpha )
    Arg.facealpha = 0.1;
end
if ~isfield( Arg, 'linewidth' )
    Arg.linewidth = get( 0, 'DefaultAxesLineWidth' );
end

n = size( y, Arg.dim );

m = mean( y, Arg.dim );
sd = std( y, 0, Arg.dim );
sem = sd / sqrt( n );

upperLim = squeeze( m + sem * icdf( 't', 0.975, n - 1 ) );
lowerLim = squeeze( m - sem * icdf( 't', 0.975, n - 1 ) );

h = fill( [ x fliplr( x ) ], [upperLim fliplr( lowerLim ) ], Arg.color, 'EdgeColor', 'none' );
set( h, 'FaceAlpha', Arg.facealpha )
set( get( get( h, 'Annotation' ), 'LegendInformation'), 'IconDisplayStyle', 'off' );
hold all
h( 2 ) = plot( x, m, 'Color', Arg.color, 'LineWidth', Arg.linewidth );

end

