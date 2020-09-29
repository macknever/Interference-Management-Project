

X1 = [1;-1;1;1];
X2 = [1;1;1;-1];
g1 = 0.1;
g2 = 0.1;
SNR1 = 10;
SNR2 = 10;

[Y1,Y2] = transmit_sym(X1,g1,X2,g2,SNR1,SNR2)

%% the function simulates two sequences go through a symmetric noisy channel
%% So what we got after this function is Y1 = X1 + g2*X2 + Z1 and Y2 = X2 + g1*X1 + Z2
%% Prof. Wang said both encoder and decoder know the g which is channel gain.
%% Z1 and Z2 are controled by SNR1 and SNR2, level of the noise.

function [Y1,Y2] = transmit_sym(X1,g1,X2,g2,SNR1,SNR2)
    Y1 = awgn(X1+g2*X2,SNR1);
    Y2 = awgn(X2+g1*X1,SNR2);
end
