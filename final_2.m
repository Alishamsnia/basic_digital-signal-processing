%%
% *Part a*

clc
close all
[n,Wn,beta,ftype] = kaiserord([0.25 7000],[1 0],[0.005 0.005],77000);
b = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
freqz(b,1,1024,10000)
sig = audioread('sig 1.WAV');
plot(sig);
title('original signal')

%%
% *part B*

% overlap-add method:َ
h = reshape(Num,[31,1]);

% Define the block size
blockSize = 100;

% Pad the input signal with zeros to ensure its length is divisible by blockSize

sig(end+1:end+blockSize-mod(length(sig),blockSize)) = 0;

% Initialize the output vector
output = zeros(length(sig)+length(h)-1,1);

% Apply the overlap-add method
for i = 1:blockSize:length(sig)
    % Extract a block from the input signal
    xBlock = sig(i:i+blockSize-1);
    
    % Filter the block using the FIR filter
    yBlock = conv(xBlock,h);
    
    % Add the filtered block to the output vector, after shifting it by an amount equal to the filter length minus one
    output(i:i+length(yBlock)-1) = output(i:i+length(yBlock)-1) + yBlock;
end

% Remove the zero-padding from the output signal
output = output(1:length(sig)+length(h)-1);
figure()
plot(output)
title('overlap-add method by L=100')
ylim([-0.3 0.3])

% overlap-save method:
h = reshape(Num,[31,1]);

% Define the block size
blockSize = 100;

% Pad the input signal with zeros to ensure its length is divisible by blockSize

sig(end+1:end+blockSize-mod(length(sig),blockSize)) = 0;

% Initialize the output vector
output = zeros(length(sig)+length(h)-1,1);

% Apply the overlap-save method
for i = 1:blockSize:length(sig)
    % Extract a block from the input signal
    xBlock = sig(i:i+blockSize-1);
    
    % Save the last M-1 samples of the block
    z = xBlock(end-length(h)+2:end);
    
    % Pad the block with zeros to ensure its length is equal to the block size plus the filter length minus one
    xBlock(end+1:end+length(h)-1) = 0;
    
    % Filter the block using the FIR filter
    yBlock = conv(xBlock,h);
    
    % Remove the first M-1 samples of the filtered block
    yBlock = yBlock(length(h):end);
    
    % Add the filtered block to the output vector
    output(i:i+blockSize-1) = yBlock(1:blockSize);
    
    % Use the last M-1 samples of the block as the first M-1 samples of the next block
    sig(i:i+length(z)-1) = z;
end

% Remove the zero-padding from the output signal
output = output(1:length(sig)+length(h)-1);
figure()
plot(output)
title('overlap-save method by L=100')
ylim([-0.3 0.3])

%%
% *Part d*

dct2_sig_fourier = fdct2(sig);
figure
plot(dct2_sig_fourier);
title('DCT type 2 with fourier');


dct1_fourier_sig = dct1_fourier(sig);
figure
plot(dct1_fourier_sig);
title('DCT type 1 with fourier');

dct2_sig = dct2(sig);
figure
plot(dct2_sig);
title('DCT type 2');


dct1_sig = dct1(sig);
figure
plot(dct1_sig);
title('DCT type 1');


% DCT type 1 normal method function
function dctCoefficients = dct1(inputSignal)
    N = length(inputSignal);
    dctCoefficients = zeros(1, N);
    
    for k = 1:N
        sum = 0;
        for n = 1:N
            sum = sum + inputSignal(n) * cos((pi/N) * (n - 0.5) * (k - 1));
        end
        if k == 1
            alpha = 1 / sqrt(N);
        else
            alpha = sqrt(2/N);
        end
        dctCoefficients(k) = alpha * sum;
    end
end

% DCT type 1 fft method function
function dctCoefficients = dct1_fourier(inputSignal)
    N = length(inputSignal);
    
    % Pad the input signal with zeros to make its length a power of 2
    paddedSignal = [inputSignal', zeros(1, 2^nextpow2(N) - N)];
    
    % Compute the DFT of the padded signal
    dftCoefficients = fft(paddedSignal);
    
    % Compute the DCT coefficients using the Fourier method
    dctCoefficients = real(dftCoefficients(1:N));
    
    % Normalize the DCT coefficients
    normalizationFactor = sqrt(2/N);
    dctCoefficients = normalizationFactor * dctCoefficients;
end

% DCT type 2 normal method function
function dctCoefficients = dct2(inputMatrix)
    [M, N] = size(inputMatrix);
    dctCoefficients = zeros(M, N);
    
    for u = 1:M
        for v = 1:N
            sum = 0;
            for x = 1:M
                for y = 1:N
                    sum = sum + inputMatrix(x, y) * cos((pi/M) * (x - 0.5) * (u - 1)) * cos((pi/N) * (y - 0.5) * (v - 1));
                end
            end
            alphaU = 1;
            alphaV = 1;
            if u == 1
                alphaU = 1 / sqrt(M);
            else
                alphaU = sqrt(2/M);
            end
            if v == 1
                alphaV = 1 / sqrt(N);
            else
                alphaV = sqrt(2/N);
            end
            dctCoefficients(u, v) = alphaU * alphaV * sum;
        end
    end
end

% DCT type 2 fft method function
function y = fdct2(x)
% Perform 2D FDCT type 2
%
% Inputs:
%   x: 2D matrix to be transformed
%
% Output:
%   y: Transformed matrix

% Get matrix size
[m, n] = size(x);

% Create transform matrix
T = zeros(m, n);
for i = 0:m-1
    for j = 0:n-1
        if i == 0
            alpha = sqrt(1/m);
        else
            alpha = sqrt(2/m);
        end
        if j == 0
            beta = sqrt(1/n);
        else
            beta = sqrt(2/n);
        end
        T(i+1,j+1) = alpha * beta * cos((pi*i*(2*j+1))/(2*n));
    end
end

% Transform rows
y1 = T * x';

% Transform columns
y = y1 * T;

end