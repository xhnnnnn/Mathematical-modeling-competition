function [fitresult, gof] = createFit1(p, E)


[xData, yData] = prepareCurveData( p, E );

% Set up fittype and options.
ft = fittype( 'poly6' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'E vs. p', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel p
% ylabel E
% grid on


