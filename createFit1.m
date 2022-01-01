function [fitresult, gof] = createFit1(rho, p)

[xData, yData] = prepareCurveData( rho, p );

% Set up fittype and options.
ft = fittype( 'poly4' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );




