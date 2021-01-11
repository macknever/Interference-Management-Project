
#include "SWSC.h"

// using namespace std;
// using namespace itpp;

// SWSC::SWSC(int bl)
// {
//     SWSC::block_length = bl;
//     SWSC::head = ones(block_length);
// }

void SWSC::regular_encode(bvec message, LDPC_Code &C, BPSK Mod, vec & regular_codeword)
{
    int block_length = C.get_nvar() - C.get_ncheck();
    int block_num = message.length() / block_length;
    for (int block_idx = 0; block_idx < block_num; block_idx++ )
    {
        // 1,fetch i-th block from the message
        bvec sub_message = message(block_idx * block_length, (block_idx + 1) * block_length - 1);
        // 2,encode the message
        bvec sub_codeword = C.encode(sub_message);
        // 3, modulate the codeword
        vec moded_sub_codeword = Mod.modulate_bits(sub_codeword);
        // 4, concatenate the codeword
        regular_codeword = concat(regular_codeword,moded_sub_codeword);        
    }
}


void SWSC::SWSC_encode(bvec message,LDPC_Code &C, BPSK Mod, vec &SWSC_codeword)
{
    int block_length = C.get_nvar() - C.get_ncheck();
    int message_length = message.length();
    int block_num = int(message_length / block_length);
    
    int Nvar = C.get_nvar();
    int n = Nvar / 2;
    // build a head vector u. as a known first block of first layer
    vec u = zeros(n);
    // vec tail = ones(n);    
    
    bvec sub_message,codeword;
    vec moded_codeword,v;
    for (int i = 0; i < block_num; i ++)
    {
        sub_message = message(i * block_length , (i + 1) * block_length - 1);
        // LDPC encode
        
        C.encode(sub_message,codeword);
        // BPSK modulation
        moded_codeword = Mod.modulate_bits(codeword);
        // cout << "i:" << i << "sub_message length:" << sub_message.length() 
        // << "moded_codeword length: " << moded_codeword.length() << endl;
        

        // superposition calculation

        // v is the vector in the second layer
        // superposition = u + 2v
        v = moded_codeword(0,n-1);
        // cout << "v length:" <<v.length() << endl;
        // cout << v(0,100) << endl;

        SWSC_codeword = concat(SWSC_codeword, 2*u + v);
           

        

        u = moded_codeword(n,2*n-1);
        
    } 

    SWSC_codeword = concat(SWSC_codeword, u);
    // cout<< "test:" << SWSC_codeword(2000,2010) << endl; 
    
    
    
    // cout << "block_length:" << block_length << endl;
    // cout << "message length:" << message.length() << endl;
    // cout << "SWSC_codeword length: " << SWSC_codeword.length() << endl;
    
}
/*
length of codeword:  n * (block_num + 1) 
n: Nvar / 2
*/

bvec SWSC::SWSC_decode_p2p(vec codeword, LDPC_Code &C, BPSK Mod, double N0)
{
    int block_length = C.get_nvar() - C.get_ncheck();    
    int Nvar = C.get_nvar();    
    int n = Nvar / 2;
    int block_num = codeword.length() / n - 1;

    // for loop
    vec sub_codeword_1,sub_codeword_2,sub_codeword;
    
    // this u should be the same as in SWSC_encode
    vec u = zeros(n);
    vec softbits;
    QLLRvec llr;
    bvec message,sub_message;

    bvec tmp_sub_codeword;
    vec tmp_sub_moded_codeword;
    C.set_exit_conditions(250);
    for (int i = 0;i < block_num; i += 1)
    {
        // first n-bits of the codeword
        sub_codeword_1 = codeword(i * n, (i + 1) * n - 1);
        // second n-bits of the codeword
        sub_codeword_2 = codeword((i + 1) * n, (i + 2) * n - 1);

        // concatenate 1st n bits and 2nd n bits to the codeword we r going to decode
        sub_codeword = concat(sub_codeword_1 - 2 * u,sub_codeword_2);

        softbits = Mod.demodulate_soft_bits(sub_codeword,N0);

        C.bp_decode(C.get_llrcalc().to_qllr(softbits),llr);


        sub_message = llr(0,block_length - 1) < 0;
        // sub_message = C.decode(softbits);
        // sub_message = C.decode(softbits);
        message = concat(message,sub_message);
        C.encode(sub_message,tmp_sub_codeword);
        tmp_sub_moded_codeword = Mod.modulate_bits(tmp_sub_codeword);
        u = tmp_sub_moded_codeword(n,2*n-1);
    }
    return message;
    
}

bvec SWSC::SWSC_decode_IAN_u_x2_v(AWGN_Channel awgn, vec codeword,LDPC_Code &C1, LDPC_Code &C2,BPSK Mod, double SNR1_dB,double SNR2_dB,double INR_dB)

{
    // parameters of message
    // block length, block numbers N,n = N/2
    int block_length1 = C1.get_nvar() - C1.get_ncheck();
    int block_length2 = C2.get_nvar() - C2.get_ncheck();

    int Nvar = C1.get_nvar();
    int n = Nvar / 2;

    int block_num = codeword.length() / n - 1;

    // parameter of signal
    // power of each signal
    // coefficient of message bits
    double p1 = pow(10,SNR1_dB/10); // power of 1st user
    double p2 = pow(10,SNR2_dB/10); // power of 2nd user    

    double INR = pow(10,INR_dB/10);

    double gain1 = sqrt(INR/p1);    // channel gain calculated from INR
    double gain2 = sqrt(INR/p2);
 
    double alpha = sqrt(5)/5;      // coefficient of signal 1
    double bet = 1;

    // test
    
    
    vec sub_codeword_1, sub_codeword_2, sub_codeword, x2;
    bvec sub_encoded_message_2;

    vec soft_x2,soft_code;
    vec sub_moded_encoded_message_2,sub_REmoded_codeword;
    vec u = zeros(n);
    vec v;

    bvec sub_message, sub_message_2,sub_REencoded_message;

    QLLRvec llr_x2,llr_code;

    bvec decoded_codeword;    


    // decoding procedure
    C1.set_exit_conditions(250);
    C2.set_exit_conditions(250);
    
    for (int i = 0; i < block_num; i++)
    {
        // Step 1
        // fetch i-th block from the codeword
        // i-th block index from i*n to i*n + n -1
        sub_codeword_1 = codeword(i * n, (i + 1) * n -1);
        sub_codeword_2 = codeword((i + 1) * n , (i + 2) * n-1);

        // Step 2 ,clean the pre_known sequence u
        // this sub_codeword is sum of 2u + v + x
        // since we decode x2(interference) first
        x2 = sub_codeword_1 - 2 * alpha * u;
 
        // use C2 to decode x2,  

        // Step 3, demodulate first n bits
        
        double N0 = awgn.get_noise();
        soft_x2 = Mod.demodulate_soft_bits(x2,N0);

        // Step 4, decode
            //soft decode
        C2.bp_decode(C2.get_llrcalc().to_qllr(soft_x2),llr_x2);
        sub_message_2 = llr_x2(0,block_length2-1) < 0;
            // hard decode
        //sub_message_2 = C2.decode(soft_x2);

        
        // Step 5, encode back to clean v
        C2.encode(sub_message_2,sub_encoded_message_2);
        // Step 6, modulate back
        sub_moded_encoded_message_2 = Mod.modulate_bits(sub_encoded_message_2);

        // Step 7, clean v
        v = x2 - bet * gain2 * sub_moded_encoded_message_2;
        v = 100 * v;
        // step 8, concatenate 
        // concatenate 2n together as a cleaned codeword for x1
        sub_codeword = concat(v,sub_codeword_2);

        // use C1 decode code

        // Step 9, demodulate N bits
        soft_code = Mod.demodulate_soft_bits(sub_codeword,N0);
        // Step 10,  decode
            // soft decode
        
        C1.bp_decode(C1.get_llrcalc().to_qllr(soft_code),llr_code);
        sub_message = llr_code(0,block_length1-1) < 0;
            // hard decode
        // sub_message = C1.decode(soft_code);

        
        // restore the sub_message to the result
        decoded_codeword = concat(decoded_codeword,sub_message);
        // caluculate i-th u for next block
        C1.encode(sub_message,sub_REencoded_message);
        sub_REmoded_codeword = Mod.modulate_bits(sub_REencoded_message);

        u = sub_REmoded_codeword(n,2*n-1);        

    }
    return decoded_codeword;


}

void SWSC::SWSC_decode_IAN_u_x2_v(AWGN_Channel awgn,vec codeword,LDPC_Code &C1, LDPC_Code &C2,BPSK Mod, double SNR1_dB,double SNR2_dB, double INR_dB,bvec &decoded_SWSC_message, bvec &decoded_IAN_message)
{
    // parameters of message
    // block length, block numbers N,n = N/2
    int block_length1 = C1.get_nvar() - C1.get_ncheck();
    int block_length2 = C2.get_nvar() - C2.get_ncheck();

    int Nvar = C1.get_nvar();
    int n = Nvar / 2;

    int block_num = codeword.length() / n - 1;

    // parameter of signal
    // power of each signal
    // coefficient of message bits
    double p1 = pow(10,SNR1_dB/10); // power of 1st user
    double p2 = pow(10,SNR2_dB/10); // power of 2nd user  

    double INR = pow(10,INR_dB/10);

    double gain1 = sqrt(INR/p1);    // channel gain calculated from INR
    double gain2 = sqrt(INR/p2);  
 
    double alpha = sqrt(5)/5;      // coefficient of signal 1
    double bet = 1;

    // test
    
    
    vec sub_codeword_1, sub_codeword_2, sub_codeword, x2;
    bvec sub_encoded_message_2;

    vec soft_x2,soft_code;
    vec sub_moded_encoded_message_2,sub_REmoded_codeword;
    vec u = zeros(n);
    vec v;

    bvec sub_message, sub_x2,sub_REencoded_message;

    QLLRvec llr_x2,llr_code;

    // bvec decoded_codeword;    


    // decoding procedure
    C1.set_exit_conditions(250);
    C2.set_exit_conditions(250);
    
    for (int i = 0; i < block_num; i++)
    {
        // Step 1
        // fetch i-th block from the codeword
        // i-th block index from i*n to i*n + n -1
        sub_codeword_1 = codeword(i * n, (i + 1) * n -1);
        sub_codeword_2 = codeword((i + 1) * n , (i + 2) * n-1);

        // Step 2 ,clean the pre_known sequence u
        // this sub_codeword is sum of 2u + v + x
        // since we decode x2(interference) first
        x2 = sub_codeword_1 - 2 * alpha * u;
 
        // use C2 to decode x2,  

        // Step 3, demodulate first n bits
        
        double N0 = awgn.get_noise();
        soft_x2 = Mod.demodulate_soft_bits(x2,N0);

        // Step 4, decode
            //soft decode
        C2.bp_decode(C2.get_llrcalc().to_qllr(soft_x2),llr_x2);
        sub_x2 = llr_x2(0,block_length2-1) < 0;
            // hard decode
        //sub_message_2 = C2.decode(soft_x2);

        decoded_IAN_message = concat(decoded_IAN_message,sub_x2);
        
        // Step 5, encode back to clean v
        C2.encode(sub_x2,sub_encoded_message_2);
        // Step 6, modulate back
        sub_moded_encoded_message_2 = Mod.modulate_bits(sub_encoded_message_2);

        // Step 7, clean v
        v = x2 - bet * gain2 * sub_moded_encoded_message_2;
        v = 10 * v;
        // step 8, concatenate 
        // concatenate 2n together as a cleaned codeword for x1
        sub_codeword = concat(v,sub_codeword_2);

        // use C1 decode code

        // Step 9, demodulate N bits
        soft_code = Mod.demodulate_soft_bits(sub_codeword,N0);

        // Step 10,  decode
            // soft decode        
        C1.bp_decode(C1.get_llrcalc().to_qllr(soft_code),llr_code);
        sub_message = llr_code(0,block_length1-1) < 0;
            // hard decode
        // sub_message = C1.decode(soft_code);

        
        // restore the sub_message to the result
        decoded_SWSC_message = concat(decoded_SWSC_message,sub_message);
        // caluculate i-th u for next block
        C1.encode(sub_message,sub_REencoded_message);
        sub_REmoded_codeword = Mod.modulate_bits(sub_REencoded_message);

        u = sub_REmoded_codeword(n,2*n-1);        

    }
}

// test the methods above
