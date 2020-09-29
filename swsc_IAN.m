%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lawrence 2020.07.16
%%% swsc & IAN two channels simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%there are block_num many blocks to transmit the message. each block has outlen(which is the length of codeword) bits to
%transit.

outlen = 2048;                  % length of codeword in each block 
block_num = 20;                 % number of blocks
simu_num = 5;                  % number of simulations


%% parameters needed
rv1 = 0;                        % Redundancy version, 0-3
rv2 = 3;
%%% redundancy version has 4 values, different value means different
%%% redundancy type. 

modulation = 'pi/2-BPSK';       % Modulation scheme, 'pi/2-BPSK' QPSK, 16QAM, 64QAM, 256QAM
nlayers = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% parameters for channel
SNR_dB_1 = 10;                  % power(dB) of 1st user
SNR_dB_2 = 10;                  % power(dB) of 2nd user
INR_dB = 8;                    % power(dB) of interference


p1 = 10^(SNR_dB_1/10);          % power of 1st user
p2 = 10^(SNR_dB_2/10);          % power of 2nd user

INR = 10^(INR_dB/10);           % power of interference

g1 = sqrt(INR/p1);              % channel gain calculated from INR
g2 = sqrt(INR/p2);

%%% This alpha is a parameter to control the power of X1 and X2
%%% cuz when generating these sequences, we produce 0s and 1s. 
alpha = sqrt(p1)*sqrt(5)/5;     % power coefficient of X1
bet = sqrt(p2);                 % power coefficient of X2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% rate of channel code

rate_list = [0.1:0.1:0.5];      % rate of 
rate_sum_list = [0.1:0.1:0.9];
rate_list_sub = [0.1:0.01:0.2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Error_list stores block error rate (BER) of each rate

Error_list_swsc_order1 = zeros(length(rate_list),1);
Error_list_swsc_order2 = zeros(length(rate_list),1);
Error_list_swsc_order3 = zeros(length(rate_list),1);
% Error_list_INA_order1 = zeros(length(rate_list),1);
% Error_list_INA_order2 = zeros(length(rate_list),1);
% Error_list_INA_order3 = zeros(length(rate_list),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% this is part 1
%%% SWSC
%%% 2 users, first using Sliding window modulation
%%% second user acts like interference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for loop for every rate in the rate_list
for column = 1:length(rate_list)
    for s = 1:simu_num
        %%% symmetric rate
        rate_1 = rate_list(column); % rate of channel code of user 1
        rate_2 = rate_1+0.05;            % rate of channel code of user 2

        %%% info_len is the length of original sequence
        %%% rate = info_len/outlen
        info_len = round2even(outlen*rate_1);

        %%% since the outlen == 2n 
        n = outlen/2;

        %%% for user2
        %%% n is the codeword length
        %%% info_len_2 is the length of original sequence of interference
        %%% transmission

        info_len_2 = round2even(n*rate_2);

        %% there are two sequence. m1 should go through swsc.m2 is regularly transmitted
        %%% m1 and m2 are both randomly created 0s and 1s. length is info_len
        %%% and info_len_2. type is double
        m1 = randi([0 1],info_len,block_num,'double');         
        m2 = randi([0 1],info_len_2,block_num,'double');


        %% parameters outlen is the length of original message. 
        %%% https://www.mathworks.com/help/5g/ref/nrdlschinfo.html?s_tid=srchtitle
        cbsInfo_1 = nrDLSCHInfo(info_len,rate_1);
        cbsInfo_2 = nrDLSCHInfo(info_len_2,rate_2);

        %% enc is for encoded message . has two levels shape like n*(blk+1) by 2 === 1024*21 by 2 .
        enc = swsc_whole_encoding(m1,outlen,cbsInfo_1,rv1);
        %%% superposition calculation x = 2u+v
        X1 = FourPam_cal(enc);    
        %%% vectorize is encode usr2's sequence into codeword length is n*(blk+1)
        X2 = vectorize(m2,n,cbsInfo_2,rv2);

        %% enpower X1 and X2
        %X1 = alpha * X1;
        %X2 = bet * X2;


        %% X1 and X2 need to be transmitted through Gaussian channel
        [Y1,Y2] = transmit_sym(X1,g1,X2,g2,alpha,bet);
        %% swsc decoding 
        % decoding order is u->x2->v
        [m1_hat,m2_hat] = correction_1(Y1,outlen,info_len,rate_1,rate_2,rv1,rv2,cbsInfo_1,cbsInfo_2,g2,alpha,bet);
        Error_list_swsc_order1(column) = Error_list_swsc_order1(column)+error_R(m1,m1_hat);
        %Error_list_INA_order1(column) = Error_list_INA_order1(column)+error_R(m2,m2_hat);
        
        [m1_hat,m2_hat] = correction_2(Y1,outlen,info_len,rate_1,rate_2,rv1,rv2,cbsInfo_1,cbsInfo_2,g2,alpha,bet);
        Error_list_swsc_order2(column) = Error_list_swsc_order2(column)+error_R(m1,m1_hat);
        %Error_list_INA_order2(column) = Error_list_INA_order2(column)+error_R(m2,m2_hat);
        
        [m1_hat,m2_hat] = correction_3(Y1,outlen,info_len,rate_1,rate_2,rv1,rv2,cbsInfo_1,cbsInfo_2,g2,alpha,bet);
        Error_list_swsc_order3(column) = Error_list_swsc_order3(column)+error_R(m1,m1_hat);
       % Error_list_INA_order3(column) = Error_list_INA_order3(column)+error_R(m2,m2_hat);
        
        fprintf('The %d th rate the %d th simulation is done \n', column,s);
    end
    Error_list_swsc_order1(column) = Error_list_swsc_order1(column)/simu_num;
    Error_list_swsc_order2(column) = Error_list_swsc_order2(column)/simu_num;
    Error_list_swsc_order3(column) = Error_list_swsc_order3(column)/simu_num;
   % Error_list_INA_order1(column) = Error_list_INA_order1(column)/simu_num;
   % Error_list_INA_order2(column) = Error_list_INA_order2(column)/simu_num;
   % Error_list_INA_order3(column) = Error_list_INA_order3(column)/simu_num;
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% simulate rate combination
for column = 1:length(rate_sum_list)
    
    %%% sum rate
    %%% looking for combination for rate which sum up for the largest that
    %%% can make the BER under a certain value like 0.1.
    sum_rate = rate_sum_list(column);
    for c = 1:length(rate_list_sub)
        rate_a = rate_list_sub(c);        
        rate_b = sum_rate - rate_a;
    end
     %%symmetric rate
    %%% rate combination

end

%% round to even
function r = round2even(var)
    var = ceil(var);
    if rem(var,2) == 0
        r = var;
    else
        r = var - 1;
    end
end
%% Block error Rate
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
%% swsc_blk_encoding is to encode sequence to a 2n == outlen codeword 
%  and split it into equal length part, 1st part is v; 2nd part is u;
function [v,u] = swsc_blk_encoding(in,outlen,cbsInfo,rv)

    tbIn = nrCRCEncode(in,cbsInfo.CRC);    
    cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN);
    ldCode = nrLDPCEncode(cbsIn,cbsInfo.BGN);
    rm_ldCode = nrRateMatchLDPC(ldCode,outlen,rv,'pi/2-BPSK',1);  
    out = double(1-2*rm_ldCode);
    %out = rm_ldCode;
    v = out(1:outlen/2);
    u = out(outlen/2+1:end);
end

%%% swsc_whole_encoding is to prepare a two layer codeword to be
%%% superposition calculated
function out = swsc_whole_encoding(in,outlen,cbsInfo,rv)
    layer_1=[];
    layer_2=[];
    head = ones(outlen/2,1);        % head is the very beginning block which we should know every single bit
                                    % here i set it to all 1s
    tail = ones(outlen/2,1);        % tial is to make the length of both layer the same
    
    [~,blks]=size(in);              % blks is the number of blocks set at the start of this snippet
    for b = 1:blks
        m = in(:,b);                % m is the sequence of this blk
        % [v u] are the encoded codeword of m.
        % v is the first half
        % u is the second half
        [v,u] = swsc_blk_encoding(m,outlen,cbsInfo,rv);   
        
        if 1==blks                  % if there is only 1 blk; only happens when testing
            layer_1 = cat(1,head,u);
            layer_2 = cat(1,v,tail);
        elseif b == 1               % otherwise, from the 1st blk
            layer_1 = cat(1,layer_1,head,u);        % concatenate the knowing head with u, make layer 1, sliding window.
            layer_2 = cat(1,layer_2,v);             % v is for layer2
        elseif b == blks                            % the last blk
            layer_1 = cat(1,layer_1,u);             % u for layer1
            layer_2 = cat(1,layer_2,v,tail);        % v and tail for layer2. here tail will never be used but to make layer1 and 2 the same length
            
        else                                        % from the 1st to the last.
            layer_1 = cat(1,layer_1,u);             % u for layer 1
            layer_2 = cat(1,layer_2,v);             % v for layer 2
        end
        
    end
   %%% layer_1 and layer_2 are both column vectors
   %%% and will be arranged from left to right
   out = cat(2,layer_1,layer_2);
end

%% vectorize encode original sequence to codeword
function codeword = vectorize(m,n,cbsInfo,rv)
    [~,c] = size(m);                                % get the size of original sequence
    codeword=[];                                    % empty codeword
    for i=1:c                                       % for each column
        info = m(:,i);                              % m is the column vector of this blk
        code = encoding(info,n,cbsInfo,rv);
        codeword = cat(1,codeword,code);
    end
    tail = ones(n,1);
    codeword = cat(1,codeword,tail);                % codeword is a length of n*(blk+1) column vector
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% the following part is standard encoding and decoding 
%%% function of LDPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sequence = decoding(y,info_len,rate,cbsInfo,rv)
    code_dec = nrRateRecoverLDPC(y,info_len,rate,rv,'pi/2-BPSK',1);
    decBits2 = nrLDPCDecode(code_dec,cbsInfo.BGN,100);               % the last argument is the number of iteration
    [blk2,blkErr] = nrCodeBlockDesegmentLDPC(decBits2,cbsInfo.BGN,info_len+cbsInfo.L);
    [sequence,tbErr] = nrCRCDecode(blk2,cbsInfo.CRC);
    
end


function codeword = encoding(in,outlen,cbsInfo,rv)
    tbIn = nrCRCEncode(in,cbsInfo.CRC); 
    cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN);
    code = nrLDPCEncode(cbsIn,cbsInfo.BGN) ;   
    c_j = nrRateMatchLDPC(code,outlen,rv,'pi/2-BPSK',1);
    codeword = double(1-2*c_j);
    %u = c_j(n+1:end)   ;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function super_position = FourPam_cal(codewords)
    u = codewords(:,1);
    v = codewords(:,2);
    super_position = [];
    for i = 1:length(u)
        %%% for each bit x = 2u+v; this is so called nature mapping
        super_position = cat(1,super_position,2*u(i)+v(i));        
    end
end

function [Y1,Y2] = transmit_sym(X1,g1,X2,g2,alpha,bet)
%     Y1 = awgn(X1+g2*X2,SNR1);
%     Y2 = awgn(X2+g1*X1,SNR2);
    
      %% N is a noise vector which power is always 1dB
    len_1 = length(X1);
    N1 = zeros(len_1,1);
   
    N = awgn(N1,0);
    
     V1 = var(X1);
     M1 = mean(X1);
        
     V2 = var(X2);
     M2 = mean(X2);
     
     Vn = var(N);
     Mn = mean(N);
     
%      %%%%%% show variance and mean of X1 and X2
%      fprintf('The variance of noise is %d and mean of noise is %d \n', Vn,Mn);
%      fprintf('The variance of X1 is %d and mean of X1 is %d \n', V1,M1);
%      fprintf('The variance of X2 is %d and mean of X2 is %d \n', V2,M2);
%      %%%%%%%% plot
%      
%      plot([X1 X2 N])
%      legend('X1','X2','noise');
    
    
   
    
    %%
    %p = 1;
    %p = 3;
    Y1 = alpha*X1 + N + bet*g2*X2;
    Y2 = bet*X2 + N + alpha*g1*X1;
end

%% this function use u->x2>v decoding order 
function [m1_hat,m2_hat] = correction_1(Y,outlen,info_len,rate_1,rate_2,rv1,rv2,cbsInfo_1,cbsInfo_2,g2,alpha,bet)%y is transmitted super_position, H the the first n bits u
    %% Y_hat is to store the prediction data.
    m1_hat=[];                              % empty decoded Y1
    m2_hat=[];                              % empty decoded Y2
    l = length(Y);
    n=outlen/2;
    b = l/n;%block_num;
    u = ones(n,1,'double');
    Y = double(Y);
       
    
    for blk = 1:b-1
        %% clean Y2
        % Y_b 1st n bits of Y in this sliding window
        % Y_b_prime 2nd n bits of Y in this sliding window
        % x2_b is n bits codeword of x2 which is y-u
        % m2_hat_b message decoded from x2_b
        % x2_b_re is n bits codeword re-encoded from m2_hat_b
        
        
        %%extract next 1:n bits from 1:2n bits of Y1
        % Y_b 
        Y_b = Y(n*blk-n+1:n*blk,1);
        %%clean out u
        x2_b = Y_b - alpha*2*double(u);
        %x2_b = x2_b/g2;
        %x2_b = x2_b/bet;
        %%% get m2 first. m2_hat_b is original sequence of m2 in this blk.
        m2_hat_b = decoding(x2_b,round2even(rate_2*n),rate_2,cbsInfo_2,rv2);
        % store it in m2_hat to compare with the original sequence.
        m2_hat = cat(2,m2_hat,m2_hat_b);
        % encode m2_hat_b to get x2_b_re, means this x2 is regenerated from
        % m2_hat_b
        x2_b_re = encoding(m2_hat_b,outlen/2,cbsInfo_2,rv2);
        
        % now we clean x2 out of v+x2, got v
        v_b = Y_b - 2 * alpha * double(u);
        v_b = v_b - g2 * bet * double(x2_b_re);
        v_b = v_b;
        % to get next block for u
        Y_b_prime = Y(blk*n+1:(blk+1)*n);
        % this sliding window is v_b concatenated with everything from next
        % bloc
        sliding_window = cat(1,v_b,Y_b_prime);
        m1_hat_b = decoding(sliding_window,info_len,rate_1,cbsInfo_1,rv1);
        m1_hat = cat(2,m1_hat,m1_hat_b);
        
        [~,u] = swsc_blk_encoding(m1_hat_b,outlen,cbsInfo_1,rv1);
        
    end
end

%% this function use u->v->x2 decoding order 
function [m1_hat,m2_hat] = correction_2(Y,outlen,info_len,rate_1,rate_2,rv1,rv2,cbsInfo_1,cbsInfo_2,g2,alpha,bet)%y is transmitted super_position, H the the first n bits u
    %% Y_hat is to store the prediction data.
    m1_hat=[];                              % empty decoded Y1
    m2_hat=[];                              % empty decoded Y2
    l = length(Y);
    n=outlen/2;
    b = l/n;%block_num;
    u = ones(n,1,'double');
    Y = double(Y);
       
    
    for blk = 1:b-1
        %% clean Y2
        
        %%extract next 1st bits from 1:2n bits of Y1
        Y_b = Y(n*blk-n+1:n*blk,1);
        %%clean out u
        v_b1 = Y_b - alpha*2*double(u);
        %x2_b = x2_b/g2;
        %x2_b = x2_b/bet;
        
        %%extract next 2nd n bits from 1:2n bits of Y1
        Y_b_prime = Y(blk*n+1:(blk+1)*n);
        
        % this sliding window is v_b concatenated with everything from next
        % bloc
        sliding_window = cat(1,v_b1,Y_b_prime);
        
        % 1st.decoding v
        v_hat_b = decoding(sliding_window,info_len,rate_1,cbsInfo_1,rv1);
        % store v to m1
        m1_hat = cat(2,m1_hat,v_hat_b);
        % regenerate codeword v to clean x2
        [v,u] = swsc_blk_encoding(v_hat_b,outlen,cbsInfo_1,rv1);
        % clean v out of [alpha*v + g*bet*x2]
        x2_b = v_b1 - alpha * double(v);
        
        %%% get m2 2nd. m2_hat_b is original sequence of m2 in this blk.
        m2_hat_b = decoding(x2_b,round2even(rate_2*n),rate_2,cbsInfo_2,rv2);
        % store it in m2_hat to compare with the original sequence.
        m2_hat = cat(2,m2_hat,m2_hat_b);
        % encode m2_hat_b to get x2_b_re, means this x2 is regenerated from
        % m2_hat_b  
        
             
    end
end

%% this function use x2->u->v decoding order
function [m1_hat,m2_hat] = correction_3(Y,outlen,info_len,rate_1,rate_2,rv1,rv2,cbsInfo_1,cbsInfo_2,g2,alpha,bet)%y is transmitted super_position, H the the first n bits u
    %% Y_hat is to store the prediction data.
    m1_hat=[];                              % empty decoded Y1
    m2_hat=[];                              % empty decoded Y2
    l = length(Y);
    n=outlen/2;
    b = l/n;%block_num;
    u = ones(n,1,'double');
    Y = double(Y);
       
    
    for blk = 1:b-1
        %%% clean Y2
        
        %%extract next 1st bits from 1:2n bits of Y
        Y_b = Y(n*blk-n+1:n*blk,1);
        
        
        %%% 1st.get m2_hat_b is original sequence of m2 in this blk.
        m2_hat_b = decoding(Y_b,round2even(rate_2*n),rate_2,cbsInfo_2,rv2);
        % store it in m2_hat to compare with the original sequence.
        m2_hat = cat(2,m2_hat,m2_hat_b);
        
        %regenerate x2
        x2_b_re = encoding(m2_hat_b,n,cbsInfo_2,rv2);
        
        % clean out u and x2 to get first n bits of v,v_b1
        v_b1 = Y_b - 2*alpha*double(u)-g2*bet*double(x2_b_re);
        %v_b1 = 4*v_b1;
        %%extract next 2nd n bits from 1:2n bits of Y1
        Y_b_prime = Y(blk*n+1:(blk+1)*n);
        
        % this sliding window is v_b concatenated with everything from next
        % bloc
        sliding_window = cat(1,v_b1,Y_b_prime);
        % decoding v 
        v_hat_b = decoding(sliding_window,info_len,rate_1,cbsInfo_1,rv1);
        
        % store v to m1
        m1_hat = cat(2,m1_hat,v_hat_b);
        % regenerate codeword v to clean x2
        [~,u] = swsc_blk_encoding(v_hat_b,outlen,cbsInfo_1,rv1);
       
             
    end
end