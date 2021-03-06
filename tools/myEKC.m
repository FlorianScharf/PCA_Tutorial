function [ nFactors, lambda, refs ] = myEKC( R, N )
% Empirical Kaiser Criterion
%
% Author: Florian Scharf, florian.scharf@uni-muenster.de and Andreas Widmann, widmann@uni-leipzig.de
% Copyright (c) 2021 Florian Scharf, University of Münster and Andreas Widmann, University of Leipzig
%
% References:
%   [1] Braeken J., & van Assen, M.A.L.M. (2017). An empirical Kaiser
%       Criterion. Psychological Methods, 22(3), 450-466. DOI:
%       10.1037/met0000074
%   [2] Steiner, M.D., & Grieder, S. (2020). EFAtools: An R package with
%       fast and flexible implementations of exploratory factor analysis
%       tools. Journal of Open Source Software, 5(53), 2521. DOI:
%       10.21105/joss.02521

p = size( R, 2 );
lambda = sort( eig( R ), 'descend' );

refs = zeros( 1, p );
for i = 1:p
    refs( i ) = max( ( ( 1 + sqrt( p / N ) ) ^ 2 ) * ( p - sum( refs ) ) / ( p - i + 1 ), 1 );
end

nFactors = find( lambda <= refs, 1 ) - 1;

end
