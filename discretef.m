clc

clear all

% This program takes an optical transmission signal and Fourier transforms it, and plots the resulting intensity


% here we insert a transmission signal to be Fourier transformed, here it is a "double slit" familiar from optics

% Slits centered at -15 and 15, width 20, height 1

for i = 1:224

    x(i) = 0;

end

for i=225:245

    x(i) = 1;

end

for i = 246:254

    x(i) = 0;

end

for i = 255:275

    x(i) = 1;

end

for i = 276:500

    x(i) = 0;

end
  
% Create horizontal axis for the graphs

t=linspace(-250,249,500);

% Fourier transform the signal and plot the DFT

y1=Ft1(x);

figure(1)
  
title('Slits of width 20, height 1 centered at +/-15')
  
plot(t,abs(y1))
  
% below is the discrete Fourier transform function

function [Y] = Ft1 ( X)


% Fourier transform function for a 1-dimensional signal.

% Input:

% X: a column vector of signal with N elements (N x 1).

% Output:

% Y: the Fourier transformed signal in the form of a column vector of N complex elements.
 

N=250; 

for u = 1:2*N % takes into account the fact that we need 500 elements in the array

    RF(u)= 0; % Initialize the output arrays

    IF(u)= 0;

    Y(u) = 0;

    for x = -N:(N-1) % here x is the x-value of the original function, the variable u must also be transformed to match

        RF(u) =RF(u) + 1/(2*N)*(real(X(x+251))*cos(-pi*x*(u-251)/N)-imag(X(x+251))*sin(-pi*x*(u-251)/N)); %  Real part of the transform

        IF(u) =IF(u) + 1/(2*N)*(real(X(x+251))*sin(-pi*x*(u-251)/N)-imag(X(x+251))*cos(-pi*x*(u-251)/N)); % Imaginary part of transform

        Y(u)=RF(u)+1i*IF(u); % the final output is the sum of real and imaginary parts

    end

end

 
end


