clc

clear all

% Convolution of transmission functions
% Here we visualise the convolution theorem in optics context
 

% Create a single slit, centred at origin, width 20, height 1

 
for i = 1:241

    x(i) = 0;

end

for i = 242:260

    x(i) = 1;

end

for i = 261:500

    x(i) = 0;

end

y1 = Ft1(x); % Fourier transform single slit

u1 = convolution1(x,x); % Convolve slit with itself

y1 = Ft1(u1); % Fourier transform the convolution

% Plot the square of the Fourier transform and the Fourier transform of the convolution

figure(1)

tiledlayout(2,1)

nexttile

plot(t,abs(y1).^2,'b')

title('Squared Fourier Transform of the Single Slit')

nexttile

plot(t,abs(y2),'r')

title('Fourier Transform of the Convolution')

% Function for the Fourier tranform

function [Y] = Ft1 ( X)


% Fourier transform function for a 1-dimensional signal.

% Input:

% X: a column vector of signal with N elements (N x 1).

%

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

% Function for the convolution

function [w] = convolution1(u , v)

% Perform convolution of two 1-dimensional signals, u and v

% Input:

% u, v: column vectors of signal with N elements each (N x 1).

% Output:

% w: the convolved signal of u and v with the same size as the input vector, N x 1


N=250;

m = length(u); % number of elements in the input arrays

n = length(v);

for k = 1:(4*N-1) % the length of the full convolved signal

    s(k) = 0;     % initialize the convolved signal

    for j = max(1,k+1-n):1:min(k,m)  % the allowed logical values for the array indices

        s(k) = s(k) + 1/(2*N)*u(j)*v(k-j+1); % discrete convolution operation, with a scaling factor 1/2N corresponding to dx

    end

end

 

for i = 1:500 % crop the signal to contain the 500 central elements

    w(i) = s(i+250);

end

end


