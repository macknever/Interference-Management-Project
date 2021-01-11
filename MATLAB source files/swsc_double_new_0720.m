%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lawrence 2020.04.16
%%% Kyle edit 2020.07.18
%%% swsc with normal encoded x2, different rv versions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%there are block_num blocks to transmit the message. each block has outlen bits to
%transit.
clc;clear all
%outlen = 1000; %%% length of each original message
block_num = 100; %%%n numbers of blocks
simu_num = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters needed
rv_1 = 0;                % Redundancy version, 0-3
rv_2 = 3;
modulation = 'pi/2-BPSK';   % Modulation scheme, QPSK, 16QAM, 64QAM, 256QAM
nlayers = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters for channel
SNR_dB_1 = 10;
SNR_dB_2 = 10;
%INR_dB = 11;


p1 = 10^(SNR_dB_1/10);
p2 = 10^(SNR_dB_2/10);

% This alpha is a parameter to control the power of X1
alpha = sqrt(p1)*sqrt(5)/5; 
% bet is for X2, SNR2 is the power is X2;
bet = sqrt(p2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rate is outlen/info_len

rate_list = [0.1:0.1:0.7];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Error_list stores error rate of each rate

Error_list_1 = zeros(length(rate_list),10);
Error_list_2 = zeros(length(rate_list),10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for loop for every rate in the rate_list

for column = 1:length(rate_list)
    col_num = 0;
    for INR_dB = 9:1:12
        INR = 10^(INR_dB/10);
        g1 = sqrt(INR/p1);
        g2 = sqrt(INR/p2);
        col_num = col_num+1;
            for i = 1:simu_num
                rate_1 = rate_list(column);
                rate_2 = 0.5;    
                
                % each original length has a corresponding codeword length.
                % info_len is codeword length
                %info_len = round2even(outlen/rate_1);  %K: codeword length is n, outlen is k, message length
                info_len = 2048;
                outlen = round2even(info_len*rate_1);
                
                
                % after encode, the codeword length is 2n and in each block, the sequence length is n.
                n = info_len/2; %K: The encoder uses a binary channel code for length 2n bits to form codeword

                % there are two message. m1 should go through swsc channel.m2 goes to normal Gaussian channel
                m1 = randi([0 1],outlen,block_num,'double');
                m2 = randi([0 1],outlen/2,block_num,'double');    

                % parameters outlen is the length of original message. 
                cbsInfo_1 = nrDLSCHInfo(outlen,rate_1);
                cbsInfo_2 = nrDLSCHInfo(outlen/2,rate_2);

                % enc is for encoded message . has two levels shape like l x 2.
                enc = swsc_whole_encoding(m1,info_len,cbsInfo_1,rv_1);       %K: should enc be +/- 1?
                % X1 is codeword after 4pam.
                X1 = FourPam_cal_alpha(enc,alpha);
                
                % X2 is normal code word
                X2 = vectorize(m2,n,cbsInfo_2,rv_2);
                X2 = bet*X2;
                [Y1,Y2] = transmit_sym(X1,g1,X2,g2,alpha,bet);

                %[Y1_hat,Y2_hat] = correction(Y1,Y2,n,info_len,rate_1,rate_2,cbsInfo_1,cbsInfo_2,g1,g2,outlen,alpha,bet,rv_1,rv_2);
                   
                [Y1_hat,Y2_hat] = new_correction(Y1,n,info_len,rate_1,rate_2,cbsInfo_1,cbsInfo_2,g1,g2,outlen,alpha,bet,rv_1,rv_2);


                Error_list_1(column,col_num) = Error_list_1(column,col_num)+error_R(m1,Y1_hat);
                Error_list_2(column,col_num) = Error_list_2(column,col_num)+error_R(m2,Y2_hat);
                fprintf('The %d th rate the %d th simulation is done \n', column,i);
            end
        Error_list_1(column,col_num) = Error_list_1(column,col_num)/simu_num;
        Error_list_2(column,col_num) = Error_list_2(column,col_num)/simu_num;
        %fprintf('INR is %.3f',g1^2*SNR1);
        fprintf('INR %.3f',INR);
    end
end

function sequence = vectorize(m,n,cbsInfo,rv_ver)
    [r,c] = size(m);
    sequence=[];
    for i=1:c
        info = m(:,i);
        codeword = encoding(info,n,cbsInfo,rv_ver);
        sequence = cat(1,sequence,codeword);
    end
    suffix = ones(n,1);
    sequence = cat(1,sequence,suffix);
end



function [Y1_hat,Y2_hat] = correction(Y1,Y2,n,info_len,rate_1,rate_2,cbsInfo_1,cbsInfo_2,g1,g2,outlen,alpha,bet,rv_1,rv_2)%y is transmitted super_position, H the the first n bits u
    %% Y_hat is to store the prediction data.
    Y1_hat=[];
    Y2_hat=[];
    l = length(Y1);
    b = l/n;%block_num;
    u = ones(n,1,'double');
    Y1 = double(Y1);
    Y2 = double(Y2);    
    
    for blk = 1:b-1
        %% clean Y2
        
        %%separate n bits from Y2
        y2 = Y1(n*blk-n+1:n*blk,1);
        %%clean out u
        x2 = y2 -alpha*2*double(u);
        x2 = x2/(bet*g2);
        x2_hat = decoding(x2,outlen/2,rate_2,cbsInfo_2,rv_2);
        
        Y2_hat = cat(2,Y2_hat,x2_hat);
        x2_hat = encoding(x2_hat,info_len/2,cbsInfo_2,rv_2);
        
        v_1 = Y1(n*blk-n+1:n*blk,1)-2*alpha*double(u);
        v_1 = v_1-g2*bet*x2_hat;
        
        y_2 = Y1(blk*n+1:(blk+1)*n);
        y = cat(1,v_1,y_2);
        y_hat = decoding(y,outlen,rate_1,cbsInfo_1,rv_1);
        Y1_hat = cat(2,Y1_hat,y_hat);
        u = encoding(y_hat,info_len,cbsInfo_1,rv_1);
        u = snd(u);
    end
end


function error_Rate = error_R(m,Y_hat)
    R=0;
    [~,blk] = size(m);
    for b=1:blk
        if isequal(m(:,b),Y_hat(:,b))
            R = R + 1;
        end
    end
    error_Rate = 1-R/blk;
end
function latter = snd(sequence)
    [len,whatever] = size(sequence);
    latter = sequence(len/2+1:end);
end

function guess = decoding(y,info_len,rate,cbsInfo,rv_ver)
    code_dec = nrRateRecoverLDPC(y,info_len,rate,rv_ver,'pi/2-BPSK',1);
    decBits2 = nrLDPCDecode(code_dec,cbsInfo.BGN,25);
    [blk2,blkErr] = nrCodeBlockDesegmentLDPC(decBits2,cbsInfo.BGN,info_len+cbsInfo.L);
    [guess,tbErr] = nrCRCDecode(blk2,cbsInfo.CRC);
    
end
function u = encoding(in,info_len,cbsInfo,rv_ver)
    tbIn = nrCRCEncode(in,cbsInfo.CRC); 
    cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN);
    code = nrLDPCEncode(cbsIn,cbsInfo.BGN) ;   
    c_j = nrRateMatchLDPC(code,info_len,rv_ver,'pi/2-BPSK',1);  
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

function [u, v] = swsc_blk_encoding(in,outlen,cbsInfo,rv_ver)

    tbIn = nrCRCEncode(in,cbsInfo.CRC);
    
    cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN);
    ldCode = nrLDPCEncode(cbsIn,cbsInfo.BGN);
    rm_ldCode = nrRateMatchLDPC(ldCode,outlen,rv_ver,'pi/2-BPSK',1);  
    out = double(1-2*rm_ldCode);
    %out = rm_ldCode;
    u = out(1:outlen/2);
    v = out(outlen/2+1:end);
end
function out = swsc_whole_encoding(in,outlen,cbsInfo,rv_ver)
    layer_1=[];
    layer_2=[];
    head = ones(outlen/2,1);
    
    [~,blks]=size(in);
    for b = 1:blks    %K: for every block
        m = in(:,b);    %K: message is each column 
        [v , u] = swsc_blk_encoding(m,outlen,cbsInfo,rv_ver);
        if 1==blks
            layer_1 = cat(1,head,u);
            layer_2 = cat(1,v,head);
        elseif b == 1
            layer_1 = cat(1,layer_1,head,u);
            layer_2 = cat(1,layer_2,v);
        elseif b == blks
            layer_1 = cat(1,layer_1,u);
            layer_2 = cat(1,layer_2,v,head);           
        else
            layer_1 = cat(1,layer_1,u);
            layer_2 = cat(1,layer_2,v);
        end
        
    end
   out = cat(2,layer_1,layer_2);
end


% new correction function
function [Y1_hat,Y2_hat] = new_correction(Y1,n,info_len,rate_1,rate_2,cbsInfo_1,cbsInfo_2,g1,g2,outlen,alpha,bet,rv_1,rv_2)
    b = length(Y1)/n;
    v = ones(n,1,'double');
    Y1 = double(Y1);
    Y1_hat=[];
    Y2_hat=[];
    for block = 1:b-1
        %First decode x2 to get Y2
        x_2 = Y1(n*block-n+1:n*block,1)-alpha*double(v);
        x_2 = x_2/bet/g2;
        x_2_hat = decoding(x_2,outlen/2,rate_2,cbsInfo_2,rv_2);
        Y2_hat = cat(2,Y2_hat,x_2_hat);
        
        %Find the decoded version of x2, subtract it from Y1
        %Then get rid of v to decode u
        decoded_x2  = encoding(x_2_hat,info_len/2,cbsInfo_2,rv_2);
        u =  Y1(n*(block-1)+1:n*block,1) - g2*decoded_x2 - alpha*v;
        x_1 = cat(1,u,Y1(block*n+1:(block+1)*n));
        x1_hat = decoding(x_1,outlen,rate_1,cbsInfo_1,rv_1);
        Y1_hat = cat(2,Y1_hat,x1_hat);
        
        %encode x1 back and find v which is the first part
        v = encoding(x1_hat,info_len,cbsInfo_1,rv_1);
        [v_len,~] = size(v); 
        v = v(1:v_len/2);
    end
end





%% This function calculate the superpositon of m and v.
%% m and v should be two l x 1 vectors. l is the length of the vector.
%% the super position is -3alpha, -alpha, alpha, 3alpha.

function super_position = FourPam_cal_alpha(codewords,alpha)
    m_j1 = codewords(:,1);
    m_j2 = codewords(:,2);
    super_position = [];
    for i = 1:length(m_j1)
        super_position = cat(1,super_position,alpha*2*m_j1(i)+alpha*m_j2(i));        
    end
end

%% the function simulates two sequences go through a symmetric noisy channel
%% So what we got after this function is Y1 = X1 + g2*X2 + Z1 and Y2 = X2 + g1*X1 + Z2
%% Prof. Wang said both encoder and decoder know the g which is channel gain.
%% Z1 and Z2 are controled by SNR1 and SNR2, level of the noise.


function [Y1,Y2] = transmit_sym(X1,g1,X2,g2,alpha,bet)
%     Y1 = awgn(X1+g2*X2,SNR1);
%     Y2 = awgn(X2+g1*X1,SNR2);
    
      %% N is a noise vector which power is always 1dB
    len_1 = length(X1);
    N = zeros(len_1,1);
    N = awgn(N,1,'measured','dB');
   % N = awgn(N,0);
    
    % N = awgn(N,1);
    
    %%
    %p = 1;
    %p = 3;
    Y1 = X1 + N + g2*X2;
    Y2 = X2 + N + g1*X1;
end