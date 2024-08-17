
function []=TransferMatrix3()

%TransferMatrix3 calculates transmission, reflection and absorption of a multilayer of planar homogenous films
%inputs:
%angle of incidence: thetai
%wavelength of incident light: lambda
%thicknesses of the layers: h
%refractive index of the layers: n (may be absorbing and dispersive, but this may require a subfunction)
%polarization: s or p (need one calculation for each for unpolarized light)

clear all
close all

thetai=50:0.01:90; %angle of incidence (degrees)
lambda=628; %vacuum wavelength (nm)
h=[NaN,42,NaN,NaN]; %film thicknesses in nm, equal in length to n, start and end with NaN
pol=1; %polarization, 1 for p and 0 for s
n=[1.6,sqrt(SR(lambda)),1.48,1.38]; %refractive index data, NaN for frequency dependence
ds=0:15:90; %film thicknesses
for b=1:length(ds)
h(3)=ds(b);
for a=1:length(thetai)
[FR(a,b),FT(a,b),FA(a,b)]=Fresnel(lambda,thetai(a),h,n,pol);
end
disp([num2str(b/length(ds)*100) '% done...'])
end
ax = axes;
ax.ColorOrder = [1 0 0; 0 0 1];
%plot results:
figure
hold on
%plot(thetai,FA)
%plot(thetai,FR,'r')
%C = {'g','b','c','m','y',[.5 .6 .7],[.8 .2 .6]};
%a=1:7;
%line_color = ['g', 'b','c','m'];
plot(thetai,FR,'LineWidth',2.5);
line_color = ["r" "b"  "k" "g" "m" "y" "c"];
colororder(line_color)
vmax = max(FR);
vmin = min(thetai);
%Find the half max value.
halfMax = max(FR) / 2;
% Find where the data first drops below half the max.
index1 = find(FR >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(FR >= halfMax, 1, 'last');
%fwhm = index2-index1 + 1 % FWHM in indexes
% if you have an x vector
fwhmFR= FR(index2) - FR(index1);
grid on;
%yticks([0, 10 20 30 40 50 60 70 80 90 100]);
yticks([0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]);
%legend('(Al_{2}O_{3}/Si_{3}N_{4})^{8}');
%'(Al_{2}O_{3}/Si_{3}N_{4})^{8}',
%xlim([])
set(gca,'XTicklabel',[800  900 1000 1100 1200],'FontSize',14)
xlabel(' \bf Wavelength (nm)','FontSize',16)
ylabel(' \bf Reflectivity','FontSize',14)
legend('N=4', 'N=6', 'N=8', 'N=10', 'N=12', 'N=14', 'N=16')
legend('boxoff')

%set(gca,'XTicklabel',[600  700   800  900 1000],'FontSize',14)
box on;
%at last 

end

function [FR,FT,FA]=Fresnel(lambda,thetai,h,n,pol)

theta(1)=thetai*pi/180;
for a=1:length(n)-1
theta(a+1)=real(asin(n(a)/n(a+1)*sin(theta(a))))-1i*abs(imag(asin(n(a)/n(a+1)*sin(theta(a)))));
end

%Fresnel coefficients:
if pol==0 %formulas for s polarization
for a=1:length(n)-1
Fr(a)=(n(a)*cos(theta(a))-n(a+1)*cos(theta(a+1)))/(n(a)*cos(theta(a))+n(a+1)*cos(theta(a+1)));
Ft(a)=2*n(a)*cos(theta(a))/(n(a)*cos(theta(a))+n(a+1)*cos(theta(a+1)));
end
elseif pol==1 %formulas for p polarization
for a=1:length(n)-1
Fr(a)=(n(a)*cos(theta(a+1))-n(a+1)*cos(theta(a)))/(n(a)*cos(theta(a+1))+n(a+1)*cos(theta(a)));
Ft(a)=2*n(a)*cos(theta(a))/(n(a)*cos(theta(a+1))+n(a+1)*cos(theta(a)));
end
end

%phase shift factors:
for a=1:length(n)-2
delta(a)=2*pi*h(a+1)/lambda*n(a+1)*cos(theta(a+1));
end

%build up transfer matrix:
M=[1,0;0,1]; %start with unity matrix
for a=1:length(n)-2
M=M*1/Ft(a)*[1,Fr(a);Fr(a),1]*[exp(-1i*delta(a)),0;0,exp(1i*delta(a))];
end
M=M*1/Ft(length(n)-1)*[1,Fr(length(n)-1);Fr(length(n)-1),1];

%total Fresnel coefficients:
Frtot=M(2,1)/M(1,1);
Fttot=1/M(1,1);

%special case of single interface:
if length(n)==2
Frtot=Fr(1);
Fttot=Ft(1);
end

%total Fresnel coefficients in intensity:
FR=(abs(Frtot))^2;
FT=(abs(Fttot))^2*real(n(length(n))*cos(theta(length(n))))/real(n(1)*cos(theta(1)));
FA=1-FR-FT;

end

function epsilon=SR(lambda)



%other parameters, worse fit to J&C but seems more accurate often:
epsiloninf=1.53;
lambdap=133;
gammap=17000;
A1=0.96;
lambda1=458;
phi1=-pi/4;
gamma1=2300;
A2=1.36;
lambda2=331;
phi2=-pi/4;
gamma2=940;

for a=1:length(lambda)
epsilon(a)=epsiloninf-1/(lambdap^2*(1/lambda(a)^2+1i/(gammap*lambda(a))))...
+A1/lambda1*(exp(phi1*1i)/(1/lambda1-1/lambda(a)-1i/gamma1)+exp(-phi1*1i)/(1/lambda1+1/lambda(a)+1i/gamma1))...
+A2/lambda2*(exp(phi2*1i)/(1/lambda2-1/lambda(a)-1i/gamma2)+exp(-phi2*1i)/(1/lambda2+1/lambda(a)+1i/gamma2));
end

end