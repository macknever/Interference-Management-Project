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
    
cbsInfo = nrDLSCHInfo(outlen,rate);

%m = randi([0 1],info_len,block_num,'double');

rate_list = [0.1:0.1:0.9];
%rate_list = [1/2]
Error_list = zeros(length(rate_list),1);

snr = 10;

u = vectorize(m2,info_len,rate,n,cbsInfo);
u

%% 1st, this function can make a r x c matrix to a rc x 1 vector
%% 2nd, the function can encode the vector to a codeword

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

function r = round2even(var)
    var = ceil(var);
    if rem(var,2) == 0
        r = var;
    else
        r = var - 1;
    end
end