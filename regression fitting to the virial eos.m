clear
clc
%the purpose of this code is to make use of regression to model
%experimental data for the fitting of the virial equation of state
% entering of experimental data 
t = [ 273.15,298.15,303.15,323.15,348.15,373.15,398.15,423.15,448.15,473.15,498.15,523.15,548.15,573.15,598.15,623.15];
B = [ -222.21,-185.8,-179.4,-156.7,-133.0,-113.6,-97.3,-83.6,-71.7,-61.5,-52.4,-44.5,-37.3,-30.9,-25.0,-19.6];
C = [10360,10600,10400,9650,8660,7720,6960,6260,5680,5290,4840,4500,4130,3860,3540,3270];
crittemp = 305.3;% critical temperature given “Introduction to Chemical Engineering Thermodynamics”, J. M. Smith, H. C. Van Ness & M. M. Abbott
redutemp = t./crittemp;% reduced temperature calculated which will be used for the model

% calculating the B values using the model. differentiate the model
% equation partially with respect to each unknown.

mat1 = ones(1,16); %values for a0 is a vector of ones

a1calc = redutemp.^-1.6; %calculate values for a1
r1=[]; %initialize a vector to store all a1 values
for n=1:1:16; %initialize a for loop to calculate all the values of a1
    r1(n)=a1calc(n);
end

a2calc = redutemp.^-4.2; %calculate values for a2
r2=[]; %initialize a vector to store all a2 values
for n = 1:1:16 %initialize a for loop to calculate all the values of a2
    r2(n)=a2calc(n);
end

mat2 = [mat1' r1' r2']; %a matrix containing all initial coefficients
mat3 = mat2\B'; %use the backslash operator to solve for the final coefficients
a0 = mat3(1);
a1 = mat3(2);
a2 = mat3(3);
calcB = a0+a1.*r1+a2.*r2; %calculated B values using the model

% calculating the C values using the model. differentiate the model
% equation partially with respect to each unknown.

mat4 = ones(1,16); %values for b0 is a vector of ones

b1calc = redutemp.^-2.8; %calculate values for b1
r3=[]; %initialize a vector to store all b1 values
for n=1:1:16; %initialize a for loop to calculate all the values of b1
    r3(n)=b1calc(n);
end

b2calc = redutemp.^-3;
r4=[]; %initialize a vector to store all b2 values
for n=1:1:16; %initialize a for loop to calculate all the values of b1
    r4(n)=b2calc(n);
end

b3calc = redutemp.^-6;
r5=[]; %initialize a vector to store all b3 values
for n=1:1:16; %initialize a for loop to calculate all the values of b3
    r5(n)=b3calc(n);
end

b4calc = redutemp.^-10.5;
r6=[]; %initialize a vector to store all b4 values
for n=1:1:16; %initialize a for loop to calculate all the values of b4
    r6(n)=b4calc(n);
end
mat5 = [mat4' r3' r4' r5' r6']; %producing a matrix with the initial coefficients
mat6 = mat5\C'; 
b0 = mat6(1);
b1 = mat6(2);
b2 = mat6(3);
b3 = mat6(4);
b4 = mat6(5);

calcC = b0+b1.*r3+b2.*r4+b3.*r5+b4.*r6 ; %calculated C values using the model

tempNEW = 340:1380; %new temperature values
redutempNEW=tempNEW./crittemp; %corresponding new reduced temperature values

%calculate B and C values for the new temperature range:

calcBnew = a0+a1.*(redutempNEW).^-1.6+a2.*(redutempNEW).^-4.2; 
calcCnew = b0+b1.*(redutempNEW).^-2.8+b2.*(redutempNEW).^-3+b3.*(redutempNEW).^-6+b4.*(redutempNEW).^-10.5; 

%plot graph of the experimental and model predictions for the second and
%third virial coefficients as functions of reduced temperature
figure(1);
plot(redutemp,B,'sq');
hold on
plot(redutempNEW,calcBnew);
xlabel('Reduced Temperature');
ylabel('B values');
title('Graph of the experimental and model predictions for the second virial coefficient');
legend('Experimental Values','Model Values');

figure(2);
plot(redutemp,C,'k*');
hold on
plot(redutempNEW,calcCnew,'r');
xlabel('Reduced Temperature');
ylabel('C values');
title('Graph of the experimental and model predictions for the third virial coefficient');
legend('Experimental Values','Model Values');

%calculations for pressure
Bderiv=-1.6*a1*crittemp^1.6.*tempNEW.^-2.6-4.2*a2*crittemp^4.2.*tempNEW.^-5.2;
Cderiv=-2.8*b1*crittemp^2.8.*tempNEW.^-3.8-3*b2*crittemp^3.*tempNEW.^-4-6*b3*crittemp^6.*tempNEW.^-7-10.5*b4*crittemp^10.5.*tempNEW.^-11.5;
volume=(2.*calcCnew-Cderiv.*tempNEW)./(Bderiv.*tempNEW-calcBnew); %see report for derivation of this equation
R=83.14; %universal gas constant
Z=1+calcBnew./volume+calcCnew./volume.^2; %virial equation of state
pressure=Z.*R.*tempNEW./volume;
critpressure=48.72; %critical pressure in bar
redupressure=pressure./critpressure;

%plot graph of reduced temperature vs. reduced pressure
figure(3);
plot(redupressure,redutempNEW,'k');
xlabel('Reduced Pressure');
ylabel('Reduced Temperature');
title('Inversion Curve');



