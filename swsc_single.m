%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lawrence 2020.03.20
%%% swsc p2p simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%there are 2000 blocks to transmit the message. each block has 4000 bits to
%transit.

outlen = 4000 %1~1023
block_num = 2000

rv = 0;                % Redundancy version, 0-3
modulation = 'pi/2-BPSK';   % Modulation scheme, QPSK, 16QAM, 64QAM, 256QAM
nlayers = 1; 




%m = randi([0 1],info_len,block_num,'double');

rate_list = [0.1:0.1:0.9];
%rate_list = [1/2]
Error_list = zeros(length(rate_list),1);

snr = 10;


for column = 1:length(rate_list)
    
    
    rate = rate_list(column);
    info_len = round2even(outlen*rate);
    m1 = randi([0 1],info_len,block_num,'double');
    m2 = randi([0 1],info_len,block_num,'double');
    
    n = outlen/2;
    
    cbsInfo = nrDLSCHInfo(info_len,rate);
    
    enc = swsc_whole_encoding(m1,outlen,rate,cbsInfo);
    
    X = FourPam_cal(enc);
    
    Y=awgn(X,snr);
    
    Y_hat = correction(Y,n,info_len,rate,cbsInfo);
    
    Error_list(column) = error_R(m1,Y_hat);

    disp('done')
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



function v_hat = block_correction(y,u)
    v_hat = double(y)-2*double(u);
end

function Y_hat = correction(Y,n,info_len,rate,cbsInfo)%y is transmitted super_position, H the the first n bits u
    Y_hat=[];
    l = length(Y);
    b = l/n;%block_num;
    u = ones(n,1,'double');
    Y = double(Y);
    
    for blk = 1:b-1
        v_1 = block_correction(Y(n*blk-n+1:n*blk,1),u);
        y_2 = Y(blk*n+1:(blk+1)*n);
        y = cat(1,v_1,y_2);
        y_hat = decoding(y,info_len,rate,cbsInfo);
        Y_hat = cat(2,Y_hat,y_hat);
        u = encoding(y_hat,info_len,rate,n,cbsInfo);
    end
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
    c_j = nrRateMatchLDPC(code,2*n,0,'pi/2-BPSK',1);  
    c_j = double(1-2*c_j);
    u = c_j(n+1:end)   ;
    
end


function super_position = FourPam_cal(codewords)
    m_j1 = codewords(:,1);
    m_j2 = codewords(:,2);
    super_position = [];
    for i = 1:length(m_j1)
        super_position = cat(1,super_position,2*m_j1(i)+m_j2(i));        
    end
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