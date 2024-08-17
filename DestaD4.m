% Create some data
x = [800 900 1000 1100 1200];
y1 = [1342 1383 1452 1593 2170  ];
y2 = [0.032 0.030 0.028 0.027 0.024];
% Create a plot with two y-axes
figure;
yyaxis left
plot(x, y1, '-diamond', 'LineWidth',2.5);
ylabel('\bf S(nm/RIU)','FontSize',14);
yyaxis right
plot(x, y2, '-diamond','LineWidth',2.5);
%ylabel('Wavelength depth');
ylabel('\bf QF(RIU^{-1})','FontSize',14);
grid on
%xlim([800 1200])
set(gca,'FontSize',14)
%xlabel(' \bf Wavelength (nm)','FontSize',16)
%ylabel(' \bf Reflectivity','FontSize',14)
legend('S', 'QF', 'N=8', 'N=10', 'N=12', 'N=14', 'N=16')
legend('boxoff')
% Add title and x-axis label
%title('Plot with Two Y-Axes');
xlabel('\bf Wavelength (nm)' ,'FontSize',14);
