#include <itpp/itcomm.h>
using namespace itpp;
using namespace std;

int main(int argc, char **argv){

    if (argc < 3){
        it_info("Usage:"<<argv[0]<<" Nvar_length (n-k) n ");
        return 1;
    }

    LDPC_Parity_Regular H;

    istringstream ss1(argv[1]);
    int Nvar;
    ss1 >> Nvar;
    istringstream ss2(argv[2]);
    int k;
    ss2 >> k;
    istringstream ss3(argv[3]);
    int n;
    ss3 >> n;

    RNG_randomize();
    {
        H.generate(Nvar,n-k, n,
               "rand",  // random unstructured matrix
               "500 10");        
        LDPC_Generator_Systematic G(&H);
        LDPC_Code C(&H,&G); 
        string NVAR,N,K;
        NVAR = argv[1];
        K = argv[2];
        N = argv[3];
        string filename = "LDPC_code_" + NVAR +"_" + K +"_"+ N + "_R.it";               
        C.save_code(filename);
    }
    RNG_randomize();
    {
        H.generate(Nvar,n-k, n,
               "rand",  // random unstructured matrix
               "500 10");        
        LDPC_Generator_Systematic G(&H);
        LDPC_Code C(&H,&G); 
        string NVAR,N,K;
        NVAR = argv[1];
        K = argv[2];
        N = argv[3];
        string filename = "LDPC_code_" + NVAR +"_" + K +"_"+ N + "_R_1.it";               
        C.save_code(filename);
    }

    

    return 0;
}