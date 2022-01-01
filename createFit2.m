function [fitresult, gof] = createFit2(theta, L)



[xData, yData] = prepareCurveData( theta, L );

% Set up fittype and options.
ft = fittype( 'poly4' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );



