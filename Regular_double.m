%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lawrence 2020.04.16
%%% swsc two channels simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%there are block_num blocks to transmit the message. each block has outlen bits to
%transit.

outlen = 2000; %%% length of each original message
block_num = 20; %%%n bumbers of blocks
simu_num = 10; %%% how many simulations needed 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters needed
rv1 = 3;
rv2 = 2;% Redundancy version, 0-3
modulation = 'QPSK';   % Modulation scheme, ，pi/2-BPSK，QPSK, 16QAM, 64QAM, 256QAM
nlayers = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters for channel
SNR_dB_1 = 10;
SNR_dB_2 = 10;
INR_dB = 12;




p1 = 10^(SNR_dB_1/10);
p2 = 10^(SNR_dB_2/10);

INR = 10^(INR_dB/10);

g1 = sqrt(INR/p1);
g2 = sqrt(INR/p2);

% g1 = 0.6;
% g2 = 0.6;
% This alpha is a parameter to control the power of X1
alpha = sqrt(p1);
% bet is for X2, SNR2 is the power is X2;
bet = sqrt(p2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rate is outlen/info_len

rate_list = [0.1:0.1:0.8];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Error_list stores error rate of each rate

Error_list_1 = zeros(length(rate_list),1);
Error_list_2 = zeros(length(rate_list),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for loop for every rate in the rate_list
for column = 1:length(rate_list)
    for i = 1:simu_num
        rate_1 = rate_list(column);
        rate_2 = rate_1+0.05;

        %% each original length has a corresponding codeword length.
        %% info_len is codeword length
        info_len = round2even(outlen*rate_1); 

        %% after encode, the codeword length is 2n and in each block, the sequence length is n.
        n = outlen;

        %% there are two message. m1 should go through swsc channel.m2 goes to normal Gaussian channel
        m1 = randi([0 1],info_len,block_num,'double');
        m2 = randi([0 1],info_len,block_num,'double');    

        %% parameters outlen is the length of original message. 
        cbsInfo_1 = nrDLSCHInfo(info_len,rate_1);
        cbsInfo_2 = nrDLSCHInfo(info_len,rate_2);

       %% X1 is normal coded
        X1 = vectorize(m1,n,rate_1,cbsInfo_1,rv1);



        %% X2 is normal coded 
        X2 = vectorize(m2,n,rate_1,cbsInfo_1,rv2);

        [Y1,Y2] = transmit_sym(X1,g1,X2,g2,alpha,bet);

        [Y1_hat,Y2_hat] = regular_correction(Y1,Y2,n,g1,g2,rate_1,rate_2,rv1,rv2,cbsInfo_1,cbsInfo_2,info_len);



        Error_list_1(column) = Error_list_1(column)+error_R(m1,Y1_hat);
        Error_list_2(column) = Error_list_2(column)+error_R(m2,Y2_hat);
        fprintf('The %d th rate the %d th simulation is done \n', column,i);
        
        
        
    end
    Error_list_1(column) = Error_list_1(column)/simu_num;
    Error_list_2(column) = Error_list_2(column)/simu_num;
end

%%%% Main part end %%%%%%%

%%%% Function definition %%%%%%%%
function sequence = vectorize(m,n,rate,cbsInfo,rv)
    [r,c] = size(m);
    sequence=[];
    for i=1:c
        info = m(:,i);
        codeword = encoding(info,n,rate,cbsInfo,rv);
        sequence = cat(1,sequence,codeword);
    end
    %suffix = ones(n,1);
    %sequence = cat(1,sequence,suffix);
end



function [Y1_hat,Y2_hat] = regular_correction(Y1,Y2,n,g1,g2,rate_1,rate_2,rv1,rv2,cbsInfo_1,cbsInfo_2,outlen)%
    %% Y_hat is to store the prediction data.
    Y1_hat=[];
    Y2_hat=[];
    l = length(Y1);
    b = l/n;%block_num;
    Y1 = double(Y1);
    Y2 = double(Y2);    
    
    for blk = 1:b
        %% clean Y2
        
        %%separate n bits from Y1 and Y2
        y1 = Y1(n*blk-n+1:n*blk,1);       
        y1_hat = decoding(y1,outlen,rate_1,cbsInfo_1,rv1);
        
        
        
        y2 = Y2(n*blk-n+1:n*blk,1); 
        
        %%this line try to clean y2 from y1.
        %%x1_hat = encoding(y1_hat,n,rate_1,cbsInfo_1);
        %%y2 = y2 - g1.*x1_hat;
        y2_hat = decoding(y2,outlen,rate_2,cbsInfo_2,rv2);
        
        Y1_hat = cat(2,Y1_hat,y1_hat);
        Y2_hat = cat(2,Y2_hat,y2_hat);
        
        
    end
end

function error_Rate = error_R(m,Y_hat)
    R=0;
    [l,blk] = size(m);
    for b=1:blk
        if isequal(m(:,b),Y_hat(:,b))
            R = R + 1;
        end
    end
    error_Rate = 1-R/blk;
end


function guess = decoding(y,info_len,rate,cbsInfo,rv)
    
    code_dec = nrRateRecoverLDPC(y,info_len,rate,rv,'QPSK',1);
    decBits2 = nrLDPCDecode(code_dec,cbsInfo.BGN,25);
    [blk2,blkErr] = nrCodeBlockDesegmentLDPC(decBits2,cbsInfo.BGN,info_len+cbsInfo.L);
    [guess,tbErr] = nrCRCDecode(blk2,cbsInfo.CRC);
    
end
function u = encoding(in,info_len,rate,cbsInfo,rv)
    [in1,in2] = size(in);
    %fprintf('Input vector size %d and %d \n', in1,in2);
    tbIn = nrCRCEncode(in,cbsInfo.CRC); 
    [tbIn1,tbIn2] = size(tbIn);
    %fprintf('tbIn vector size %d and %d \n', tbIn1,tbIn2);
    cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN);
    [cbsIn1,cbsIn2] = size(cbsIn);
   % fprintf('cbsIn vector size %d and %d \n', cbsIn1,cbsIn2);
    code = nrLDPCEncode(cbsIn,cbsInfo.BGN) ;   
    [code1,code2] = size(code);
    %fprintf('code vector size %d and %d \n', code1,code2);
    c_j = nrRateMatchLDPC(code,info_len,rv,'QPSK',1) ;
    [c_j1,c_j2] = size(c_j);
    %fprintf('c_j vector size %d and %d \n', c_j1,c_j2);
    u = double(1-2*c_j);
    %u = c_j(n+1:end)   ;
    
end

function r = round2even(var)
    var = ceil(var);
    if rem(var,2) == 0
        r = var;
    else
        r = var - 1;
    end
end







%% the function simulates two sequences go through a symmetric noisy channel
%% So what we got after this function is Y1 = X1 + g2*X2 + Z1 and Y2 = X2 + g1*X1 + Z2
%% Prof. Wang said both encoder and decoder know the g which is channel gain.
%% Z1 and Z2 are controled by SNR1 and SNR2, level of the noise.


function [Y1,Y2] = transmit_sym(X1,g1,X2,g2,alpha,bet)
%       Y1 = X1 + 0.999999999999*X2;
%       Y2 = X2 + 0.9*X1;
    
%     Y1 = awgn(X1,SNR1)+g2*X2;
%     Y2 = awgn(X2,SNR2)+g1*X1;
    
%     Y1 = awgn(X1,SNR1,'measured')+g2*X2;
%     Y2 = awgn(X2,SNR2,'measured')+g1*X1;

    %% N is a noise vector which power is always 1dB
    len_1 = length(X1);
    N = zeros(len_1,1);
    N = awgn(N,0,'measured');
    
    % N = awgn(N,1);
    
    %%
    %p = 1;
      Y1 = alpha*X1 + N + bet*g2*X2;
      Y2 = bet*X2 + N + alpha*g1*X1;
   
%     
%     Y1 = awgn(X1,SNR1,'measured','linear')+g2*X2;
%     Y2 = awgn(X2,SNR2,'measured','linear')+g1*X1;
end