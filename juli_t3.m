function [fitresult, gof] = juli_t3(penyouzui_t3, penyouzui_h3)


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( penyouzui_t3, penyouzui_h3 );

% Set up fittype and options.
ft = fittype( 'poly4' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );



