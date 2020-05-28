
time50 = [3.23 13.69 22.188 140.2 654.18 3768.7];
NumberNodes = [16 64 80 160 320 480];

col=hsv(2);
marker = ['o','s','d','*','v'];

%%Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( NumberNodes, time50 );

% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [18.318765494261 0.0111171197048772];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

xAxis = linspace(0,500,500);
figure(1)
ydata = time50;
xdata = NumberNodes;
p1=plot(xdata,ydata,'Color',col(1,:));
p1.LineStyle  = 'none'; p1.LineWidth = 0.5; p1.Marker = marker(1); 
p1.MarkerSize = 8; p1.MarkerFaceColor = col(1,:);
hold on
p2 = plot(xAxis,fitresult(xAxis),'Color',col(1,:));
p2.LineStyle  = '-'; p2.LineWidth = 1.5; p2.Marker = 'none'; 

% p2 = plot(xAxis,exp(xAxis),'Color',col(2,:));
% p2.LineStyle  = '-'; p2.LineWidth = 1.5; p2.Marker = 'none'; 
% 

%ylim([0,500])

set(gca,'FontSize',14); %set the axis font size
box on;
grid on;
xlabel('Number of Nodes', 'FontSize', 18)
ylabel('time [s]', 'FontSize', 18)
hold off
title('Time to run 50 runs','FontSize', 18)
