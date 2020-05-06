%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lawrence 2020.04.16
%%% swsc two channels simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%there are block_num blocks to transmit the message. each block has outlen bits to
%transit.

outlen = 4000; %%% length of each original message
block_num = 2000; %%%n bumbers of blocks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters needed
rv = 0;                % Redundancy version, 0-3
modulation = 'pi/2-BPSK';   % Modulation scheme, QPSK, 16QAM, 64QAM, 256QAM
nlayers = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters for channel
g1 = 0.6;
g2 = 0.6;
SNR1 = 10;
SNR2 = 10;
alpha = sqrt(5)/5; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rate is outlen/info_len

rate_list = [0.18:0.01:0.23];
%rate_list = [1/2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Error_list stores error rate of each rate

Error_list_1 = zeros(length(rate_list),1);
Error_list_2 = zeros(length(rate_list),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for loop for every rate in the rate_list
for column = 1:length(rate_list)
    for i = 1:1000
        rate_1 = rate_list(column);
        rate_2 = 2*rate_1;

        %% each original length has a corresponding codeword length.
        %% info_len is codeword length
        info_len = round2even(outlen/rate_1); 

        %% after encode, the codeword length is 2n and in each block, the sequence length is n.
        n = info_len/2;

        %% there are two message. m1 should go through swsc channel.m2 goes to normal Gaussian channel
        m1 = randi([0 1],outlen,block_num,'double');
        m2 = randi([0 1],outlen,block_num,'double');    

        %% parameters outlen is the length of original message. 
        cbsInfo_1 = nrDLSCHInfo(outlen,rate_1);
        cbsInfo_2 = nrDLSCHInfo(outlen,rate_2);

       %% X1 is normal coded
        X1 = vectorize(m1,n,rate_2,cbsInfo_2);



        %% X2 is normal coded 
        X2 = vectorize(m2,n,rate_2,cbsInfo_2);

        [Y1,Y2] = transmit_sym(X1,g1,X2,g2,SNR1,SNR2);

        [Y1_hat,Y2_hat] = regular_correction(Y1,Y2,n,g1,g2,rate_2,cbsInfo_2,outlen);



        Error_list_1(column) = Error_list_1(column)+error_R(m1,Y1_hat);
        Error_list_2(column) = Error_list_2(column)+error_R(m2,Y2_hat);
        fprintf('The %d th rate the %d th simulation is done \n', column,i);
        
        
        
    end
    Error_list_1(column) = Error_list_1(column)/1000;
    Error_list_2(column) = Error_list_2(column)/1000;
end

function sequence = vectorize(m,n,rate,cbsInfo)
    [r,c] = size(m);
    sequence=[];
    for i=1:c
        info = m(:,i);
        codeword = encoding(info,n,rate,cbsInfo);
        sequence = cat(1,sequence,codeword);
    end
    suffix = ones(n,1);
    sequence = cat(1,sequence,suffix);
end

function v_hat = block_correction(y,u,alpha)
    v_hat = double(y)-2*alpha*double(u);
end

function [Y1_hat,Y2_hat] = regular_correction(Y1,Y2,n,g1,g2,rate_2,cbsInfo_2,outlen)%
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
        y1_hat = decoding(y1,outlen,rate_2,cbsInfo_2);
        x1_hat = encoding(y1_hat,n,rate_2,cbsInfo_2);
        
        
        y2 = Y2(n*blk-n+1:n*blk,1);  
        
        y2 = y2 - g1.*x1_hat;
        y2_hat = decoding(y2,outlen,rate_2,cbsInfo_2);
        
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
function latter = snd(sequence)
    [len,whatever] = size(sequence);
    latter = sequence(len/2+1:end);
end

function guess = decoding(y,info_len,rate,cbsInfo)
    code_dec = nrRateRecoverLDPC(y,info_len,rate,0,'pi/2-BPSK',1);
    decBits2 = nrLDPCDecode(code_dec,cbsInfo.BGN,25);
    [blk2,blkErr] = nrCodeBlockDesegmentLDPC(decBits2,cbsInfo.BGN,info_len+cbsInfo.L);
    [guess,tbErr] = nrCRCDecode(blk2,cbsInfo.CRC);
    
end
function u = encoding(in,info_len,rate,cbsInfo)
    tbIn = nrCRCEncode(in,cbsInfo.CRC); 
    cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN);
    code = nrLDPCEncode(cbsIn,cbsInfo.BGN) ;   
    c_j = nrRateMatchLDPC(code,info_len,0,'pi/2-BPSK',1);  
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

function [u v] = swsc_blk_encoding(in,outlen,rate,cbsInfo)

    tbIn = nrCRCEncode(in,cbsInfo.CRC); 
    cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN);
    ldCode = nrLDPCEncode(cbsIn,cbsInfo.BGN);
    rm_ldCode = nrRateMatchLDPC(ldCode,outlen,0,'pi/2-BPSK',1);  
    out = double(1-2*rm_ldCode);
    u = out(1:outlen/2);
    v = out(outlen/2+1:end);
end
function out = swsc_whole_encoding(in,outlen,rate,cbsInfo)
    layer_1=[];
    layer_2=[];
    head = ones(outlen/2,1);
    
    [~,blks]=size(in);
    for b = 1:blks
        m = in(:,b);
        [v u] = swsc_blk_encoding(m,outlen,rate,cbsInfo);
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


function [Y1,Y2] = transmit_sym(X1,g1,X2,g2,SNR1,SNR2)
    Y1 = awgn(X1+g2*X2,SNR1);
    Y2 = awgn(X2+g1*X1,SNR2);
end