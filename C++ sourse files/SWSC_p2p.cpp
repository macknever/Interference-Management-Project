#include <itpp/itcomm.h>
#include <sstream>
#include "SWSC.h"

using namespace itpp;
using namespace std;




int main(int argc, char **argv)
{
    // check arguments' validity
    if (argc < 2){
        it_info("Usage: "<<argv[0]<<" codec_file.it [SNR_dB]  [block_num]");
        return 1;
    }

    LDPC_Generator_Systematic G;   
    LDPC_Code C(argv[1],&G);  // SWSC LDPC code
    

    C.set_exit_conditions(250);
    
    BPSK Mod;
    SWSC swsc;

    
    int block_num;    

    {
        istringstream ss(argv[3]);        
        ss >> block_num;
    }
    
    double SNR_dB;

    {
        istringstream ss(argv[2]);        
        ss >> SNR_dB;
    }

    double N0 = pow(10,-SNR_dB/10);
    
    
    
    int SWSC_block_length = C.get_nvar() - C.get_ncheck();
    

    int SWSC_message_length = block_num * SWSC_block_length;
    


    bvec SWSC_message = randb(SWSC_message_length);    
    
    bvec decoded_message;
   
    vec SWSC_codeword;
    vec received_codeword;
    
    swsc.SWSC_encode(SWSC_message,C,Mod,SWSC_codeword);
    
    RNG_randomize();
    AWGN_Channel awgn(N0);
    vec p2pSWSC_out = awgn(SWSC_codeword);
    bvec p2pDecoded_message = swsc.SWSC_decode_p2p(p2pSWSC_out,C,Mod,N0); 

    cout<<"length of orginal message: " << SWSC_message.length()<< endl;
    cout<<"length of decoded message: " << p2pDecoded_message.length()<< endl;
   
    
    
    
    BLERC blerc;
    BERC berc;

    cout<<SWSC_block_length<<'/'<<SWSC_message.length()<<'/'<<p2pDecoded_message.length()<<endl;
    blerc.set_blocksize(SWSC_block_length);
    
    blerc.count(SWSC_message,p2pDecoded_message);
    
    berc.count(SWSC_message,p2pDecoded_message);

    
    
    cout << "block error rate p2p: " <<blerc.get_errorrate() << endl;
    cout << "bit error rate: " <<berc.get_errorrate() << endl;
    
    return 0;
}
