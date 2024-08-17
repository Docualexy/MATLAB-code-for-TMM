% Create some data
x = [800 900 1000 1100 1200];
y1 = [500 344.8 291.56 264.43 248.05];
y2 = [0.769 0.788 0.799 0.808 0.815];
% Create a plot with two y-axes
figure;
yyaxis left
plot(x, y1, '-diamond', 'LineWidth',2.5);
ylabel('\bf FOM(RIU)','FontSize',14);
yyaxis right
plot(x, y2, '-diamond','LineWidth',2.5);
%ylabel('Wavelength depth');
ylabel('\bf DA','FontSize',14);
grid on
%xlim([800 1200])
set(gca,'FontSize',14)
%xlabel(' \bf Wavelength (nm)','FontSize',16)
%ylabel(' \bf Reflectivity','FontSize',14)
legend('FOM', 'DA', 'N=8', 'N=10', 'N=12', 'N=14', 'N=16')
legend('boxoff')
% Add title and x-axis label
%title('Plot with Two Y-Axes');
xlabel('\bf Wavelength (nm)' ,'FontSize',14);
