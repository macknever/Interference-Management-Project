outlen = 40 % length of each original message
block_num = 4 % just number of blocks

% some parameter needed
rv = 0;                % Redundancy version, 0-3
modulation = 'pi/2-BPSK';   % Modulation scheme, QPSK, 16QAM, 64QAM, 256QAM
nlayers = 1; 


rate = 0.1
info_len = round2even(outlen/rate);
m1 = randi([0 1],outlen,block_num,'double');
m2 = randi([0 1],outlen,block_num,'double');
    
n = info_len/2;
alpha = sqrt(5)/5;    
cbsInfo = nrDLSCHInfo(outlen,rate);

%m = randi([0 1],info_len,block_num,'double');

rate_list = [0.1:0.1:0.9];
%rate_list = [1/2]
Error_list = zeros(length(rate_list),1);

snr = 10;

u = swsc_whole_encoding(m1,info_len,rate,cbsInfo);
X1 = FourPam_cal_alpha(u,alpha);
X2 = vectorize(m2,info_len,rate,n,cbsInfo);

g1 = 0.1
g2 = 0.1;
SNR1 = 10;
SNR2 = 10;

[Y1,Y2] = transmit_sym(X1,g1,X2,g2,SNR1,SNR2);
Y1
Y2


function v_hat = block_correction_(y,u,alpha,channel_gain_2,y2_hat)
    v_hat = double(y) - 2*double(u) / alpha - channel_gain_2*y2_hat;
    v_hat = v_hat/alpha;
end


function sequence = vectorize(m,info_len,rate,n,cbsInfo)
    [r,c] = size(m);
    sequence=[];
    for i=1:c
        info = m(:,i);
        codeword = encoding(info,info_len,rate,n,cbsInfo)
        sequence = cat(1,sequence,codeword);
    end
    suffix = ones(n,1);
    sequence = cat(1,sequence,suffix);
end

function guess = decoding(y,info_len,rate,cbsInfo)
    code_dec = nrRateRecoverLDPC(y,info_len,rate,0,'pi/2-BPSK',1);
    decBits2 = nrLDPCDecode(code_dec,cbsInfo.BGN,25);
    [blk2,blkErr] = nrCodeBlockDesegmentLDPC(decBits2,cbsInfo.BGN,info_len+cbsInfo.L);
    [guess,tbErr] = nrCRCDecode(blk2,cbsInfo.CRC);
    
end

function u = encoding(in,info_len,rate,n,cbsInfo)
    tbIn = nrCRCEncode(in,cbsInfo.CRC); 
    cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN);
    code = nrLDPCEncode(cbsIn,cbsInfo.BGN) ;   
    c_j = nrRateMatchLDPC(code,n,0,'pi/2-BPSK',1);  
    u = double(1-2*c_j);
    %u = c_j(n+1:end)   ;
    
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

function r = round2even(var)
    var = ceil(var);
    if rem(var,2) == 0
        r = var;
    else
        r = var - 1;
    end
end

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