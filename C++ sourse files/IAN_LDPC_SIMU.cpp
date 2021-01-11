/*This simulation is about 2 regular transmission which will interfere each other
and one always see the other interference as noise

y1 = alpha * x1 + gain2 * bet * x2 + N
y2 = bet * x2 + gain1 * alpha * x1 + N
*/




#include <itpp/itcomm.h>
#include <sstream>
#include "SWSC.h"

using namespace itpp;
using namespace std;


int main(int argc, char **argv)
{
    // check arguments' validity
    if (argc < 2){
        it_info("Usage:"<<argv[0]<<"codec_file1.it  codec_file2 [INR_dB].it");
        return 1;
    }

    double INR_dB;
    
    {
        istringstream ss(argv[3]);        
        ss >> INR_dB;
    }
   
    // build a LDPC code from existing generator and parity check matrix     
    LDPC_Generator_Systematic G1;    
    LDPC_Code C1(argv[1],&G1);
    C1.set_exit_conditions(250);

    LDPC_Generator_Systematic G2;    
    LDPC_Code C2(argv[2],&G2);
    C2.set_exit_conditions(200);

    // calculate the each code's block length
    int block_length1 = C1.get_nvar() - C1.get_ncheck();
    int block_length2 = C2.get_nvar() - C2.get_ncheck();

    // set block number
    int block_num = 30;

    // init 2 random message

    bvec message1 = randb(block_num * block_length1);
    bvec message2 = randb(block_num * block_length2);

    bvec codeword1,codeword2;

    BLERC blerc1,blerc2;
    blerc1.set_blocksize(block_length1);
    blerc2.set_blocksize(block_length2);
    // set parameters of channel and the power of signal
    

    QLLRvec llr1,llr2;              // for soft_decoding
    // modulation types
    BPSK Mod;
    
    // power of signal
    double SNR1_dB = 10;            // power(dB) of 1st user
    double SNR2_dB = 10;            // power(dB) of 2nd user

    // double INR_dB = 0;              // power(dB) of interference


    // double N0_dB = pow(10,-SNR1_dB/10);
    double N0 = 1;
    cout << "N0:" << N0 << endl;
    // double N0_dB = pow(10,-SNR1_dB/10);
    AWGN_Channel awgn(N0);
    // calculation
    double p1 = pow(10,SNR1_dB/10); // power of 1st user
    double p2 = pow(10,SNR2_dB/10); // power of 2nd user

    double INR = pow(10,INR_dB/10);

    double gain1 = sqrt(INR/p1);    // channel gain calculated from INR
    double gain2 = sqrt(INR/p2);

    //gain1 = gain2 = 0;
    
    // as the power of noise will be constant 1, so we need to set the coefficent of signals
    // cuz after encoded they are all -1s and 1s

    // alpha represents coefficient of x1
    // bet represents coefficient of x2
    // the representation above would be fixed through all simulations

    // double alpha = sqrt(p1)*sqrt(5)/5;     // 4PAM coefficient of x1 
    double alpha = sqrt(p1);                  // BPSK coefficient of x1
    double bet = sqrt(p2);                    // coefficient of x2

    // checking parameters

    cout << "p1,p2:" << p1 << "," << p2 << endl;
    cout << "INR:" << INR << endl;
    cout << "gain1,gain2:" << gain1 << "," << gain2 << endl;
    cout << "alpha,bet:" << alpha << "," << bet << endl;
    // Step 1, encoding
    // block by block
    bvec block_codeword1,block_codeword2;
    for (int k = 0;k < block_num; k++)
    {
        C1.encode(message1(k * block_length1,(k+1) * block_length1-1),block_codeword1);
        codeword1 = concat(codeword1,block_codeword1);

        C2.encode(message2(k * block_length2,(k+1) * block_length2-1),block_codeword2);
        codeword2 = concat(codeword2,block_codeword2);
    }
    

    // Step 2, modulation    
    vec modulated_codeword1 = Mod.modulate_bits(codeword1);
    vec modulated_codeword2 = Mod.modulate_bits(codeword2);

    // Step 3, randomnize
    

    // Step 4, transmission
    // y1 = alpha * x1 + gain2 * bet * x2 + N
    // y2 = bet * x2 + gain1 * alpha * x1 + N
    vec received_codeword1 = 
        awgn(alpha * modulated_codeword1 + gain2 * bet * modulated_codeword2);

    vec received_codeword2 = 
        awgn(bet * modulated_codeword2 + gain1 * alpha * modulated_codeword1);

    // Step 5, demodulation
    vec demoded_received_codeword1 = Mod.demodulate_soft_bits(received_codeword1,N0);
    vec demoded_received_codeword2 = Mod.demodulate_soft_bits(received_codeword2,N0); 

    // Step 6, decoding
    // block by block
    bvec decoded_codeword1,decoded_codeword2;
    
    for (int i = 0; i < block_num;i++)
    {
        
        // soft decode

        // C1.bp_decode(C1.get_llrcalc().to_qllr(
        //     demoded_received_codeword1(i * C1.get_nvar(),(i+1) * C1.get_nvar()-1)),llr1);
        // decoded_codeword1 = concat(decoded_codeword1,llr1(0,block_length1-1)<0);

        // C2.bp_decode(C2.get_llrcalc().to_qllr(
        //     demoded_received_codeword2(i * C2.get_nvar(),(i+1) * C2.get_nvar()-1)),llr2);
        // decoded_codeword2 = concat(decoded_codeword2,llr2(0,block_length2-1)<0);

        // hard decode
        bvec tmp_de1 = C1.decode(
            demoded_received_codeword1(i * C1.get_nvar(),(i+1) * C1.get_nvar()-1));
        decoded_codeword1 = concat(decoded_codeword1,tmp_de1);

        bvec tmp_de2 = C2.decode(
            demoded_received_codeword2(i * C2.get_nvar(),(i+1) * C2.get_nvar()-1));
        decoded_codeword2 = concat(decoded_codeword2,tmp_de2);
        

        
    }
    
    // Step 7, check block error rate
    blerc1.count(message1,decoded_codeword1);
    blerc2.count(message2,decoded_codeword2);

    cout << "block error rate of message1: " << blerc1.get_errorrate() << endl;
    cout << "block error rate of message2: " << blerc2.get_errorrate() << endl;
    

    return 0;
}