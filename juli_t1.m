function [fitresult, gof] = juli_t1(penyouzui_t1, penyouzui_h1)

[xData, yData] = prepareCurveData( penyouzui_t1, penyouzui_h1 );

% Set up fittype and options.
ft = fittype( 'poly4' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );


