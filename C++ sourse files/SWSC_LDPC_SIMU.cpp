#include <itpp/itcomm.h>
#include <sstream>
#include "SWSC.h"

using namespace itpp;
using namespace std;

vec transfer(AWGN_Channel awgn, vec codeword1,vec codeword2, double INR_dB, double SNR1_dB,double SNR2_dB)
{
    
    // double N0 = 1;
    // AWGN_Channel awgn(N0);
    // calculation
    double p1 = pow(10,SNR1_dB/10); // power of 1st user
    double p2 = pow(10,SNR2_dB/10); // power of 2nd user

    double INR = pow(10,INR_dB/10);

    double gain1 = sqrt(INR/p1);    // channel gain calculated from INR
    double gain2 = sqrt(INR/p2);
 
    double alpha = sqrt(5)/5;
    double bet = 1; 

    double N0 = awgn.get_noise();

    cout << "p1: " << p1 << endl;
    cout << "p2: " << p2 << endl;
    cout << "INA: " << INR << endl;
    cout << "gain1: " << gain1 << endl;
    cout << "gain2: " << gain2 << endl;
    cout << "N0: " << N0 << endl;

    // y = alpha * x1 + g2 * bet * x2 + N
    vec codeword = alpha * codeword1 + gain1* bet * codeword2;
    vec received_codeword = awgn(codeword);

    return codeword;
}


int main(int argc, char **argv)
{
    // check arguments' validity
    if (argc < 2){
        it_info("Usage: "<<argv[0]<<" codec_file1.it  codec_file2.it [INR_dB] [block_num]");
        return 1;
    }

    LDPC_Generator_Systematic G1,G2;    
    LDPC_Code C1(argv[1],&G1);  // SWSC LDPC code
    LDPC_Code C2(argv[2],&G2);  // IAN LDPC code

    C1.set_exit_conditions(250);
    C2.set_exit_conditions(250);
    BPSK Mod;
    SWSC swsc;

    double INR_dB;
    int block_num;
    {
        istringstream ss(argv[3]);        
        ss >> INR_dB;
    }

    {
        istringstream ss(argv[4]);        
        ss >> block_num;
    }
    // cout << "INR_dB:" << INR_dB << endl;
    cout << "block_num:" << block_num << endl;

    // parameter of signals

    double SNR1_dB = 10;            // power(dB) of 1st user
    double SNR2_dB = 10;            // power(dB) of 2nd user
    double N0 = pow(10,-SNR1_dB/10);                  // power(dB) of noise
    N0 = 1;
    double p1 = pow(10,SNR1_dB/10); // power of 1st user
    double p2 = pow(10,SNR2_dB/10); // power of 2nd user

    double INR = pow(10,INR_dB/10);

    double gain1 = sqrt(INR/p1);    // channel gain calculated from INR
    double gain2 = sqrt(INR/p2);

    
    // as the power of noise will be constant 1, so we need to set the coefficent of signals
    // cuz after encoded they are all -1s and 1s

    // alpha represents coefficient of x1
    // bet represents coefficient of x2
    // the representation above would be fixed through all simulations

    double alpha_4PAM = sqrt(p1)*sqrt(5)/5;     // 4PAM coefficient of x1 
    double alpha_BPSK = sqrt(p1);                  // BPSK coefficient of x1
    double bet_BPSK = sqrt(p2);
    
    
    int SWSC_block_length = C1.get_nvar() - C1.get_ncheck();
    int IAN_block_length = C2.get_nvar() - C2.get_ncheck();

    int SWSC_message_length = block_num * SWSC_block_length;
    int IAN_message_length = (block_num + 1) * IAN_block_length;

    RNG_randomize();
    bvec SWSC_message = randb(SWSC_message_length);
    
    RNG_randomize();
    bvec IAN_message = randb(IAN_message_length);
    
    bvec decoded_SWSC_message,decodec_IAN_message;

    vec SWSC_codeword , IAN_codeword;
    vec received_codeword;
    
    swsc.SWSC_encode(SWSC_message,C1,Mod,SWSC_codeword);
    swsc.regular_encode(IAN_message,C2,Mod,IAN_codeword);

    AWGN_Channel awgn(N0 / 2);
    vec p2pSWSC_out = awgn(SWSC_codeword);
    bvec p2pDecoded_message = swsc.SWSC_decode_p2p(p2pSWSC_out,C1,Mod,N0);

    received_codeword = transfer(awgn,SWSC_codeword,IAN_codeword,INR_dB,SNR1_dB,SNR2_dB);
    // cout << SWSC_codeword.length() << "/" << IAN_codeword.length() << endl;
    swsc.SWSC_decode_IAN_u_x2_v(awgn, received_codeword,C1,C2,Mod,SNR1_dB,SNR2_dB,INR_dB,decoded_SWSC_message,decodec_IAN_message);


    // cout<<SWSC_message.length()<<'/'<< p2pDecoded_message.length() << endl;
    
    
    BLERC blerc,blerc_p2p,blerc_IAN;
    BERC berc,berc_p2p,berc_IAN;

    blerc.set_blocksize(C1.get_nvar() - C1.get_ncheck());
    blerc_p2p.set_blocksize(C1.get_nvar() - C1.get_ncheck());
    blerc_IAN.set_blocksize(C2.get_nvar() - C2.get_ncheck());

    blerc.count(SWSC_message,decoded_SWSC_message);
    blerc_p2p.count(SWSC_message,p2pDecoded_message);
    blerc_IAN.count(IAN_message,decodec_IAN_message);

    berc.count(SWSC_message,decoded_SWSC_message);
    berc_p2p.count(SWSC_message,p2pDecoded_message);
    berc_IAN.count(IAN_message,decodec_IAN_message);

    cout << "block error rate SWSC: " <<blerc.get_errorrate() << endl;
    cout << "block error rate p2p: " <<blerc_p2p.get_errorrate() << endl;
    cout << "block error rate IAN: " <<blerc_IAN.get_errorrate() << endl;

    cout << "bit error rate SWSC: " <<berc.get_errorrate() << endl;
    cout << "bit error rate p2p: " <<berc_p2p.get_errorrate() << endl;
    cout << "bit error rate IAN: " <<berc_IAN.get_errorrate() << endl;
    
    return 0;
}
