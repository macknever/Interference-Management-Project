# README
## Introduction
Interference Management Project's source files  
There are c++ source files written by senior communication engineers and Matlab source files written by Lawrence and Kyle.   
This instruction will tell you:
- [What is this project about?](#interference-management-project)
- [About C++ source files](#cpp-source-files)
- [About Matlab source files](#matlab-source-files)
- [How to take advantage of Sockeye](#ubc-sockeye-server)
## Interference Management Project
To be continued
## CPP source files
### Compiling C++ in Linux 
These `.cpp` file should be executed under Linux. You can find an old spare laptop to install a Linux, or install VM into your Windows computer, or you can install subsystem into your Windows.  
It is not hard to install linux, there are bunch of videos on [youtube](https://www.youtube.com/?gl=CA) teaching how to install Linux. Here are instructions for how to install Itpp library after you get your Linux system.
- #### *Subsystem Ubuntu 16.04 or 20.04*
    - Intructions can be found on [IT++](http://itpp.sourceforge.net/4.3.1/installation.html) website.
    - Ubuntu 18.04 is not tested for installing or compiling.
    ```shell
    $ wget https://sourceforge.net/projects/itpp/files/itpp/4.3.1/itpp-4.3.1.tar.bz2
    $ bzip2 -cd itpp-4.3.1.tar.bz2 | tar xf -
    $ cd itpp-4.3.1
    $ mkdir build
    $ cd build
    // install cmake and g++ in case they are not installed
    $ sudo apt update
    $ sudo apt install cmake
    $ sudo apt install g++
    $ sudo apt install libitpp-dev
    ```
    Right now you can compile `.cpp`file which using functions from itpp library. To test that, you can draft a test code.
    ```shell
    $ cd $HOME/<YOUR_USER_NAME>/itpptest
    $ vim itpptest.cpp
    ```
    In this `itpptest.cpp`file type in:
    ```c++
    #include<itpp/itcomm.h>
    #include<iostream>
    
    using namespace itpp;
    
    using std::cout;
    using std::endl;
    
    int main(){
        BPSK bpsk;
        cout<<"Hello it++ library"<<endl;
        return 0;
    }
    ```
    Next, is the compiling and linking part, 
    ```shell
    $ g++ `itpp-config --cflags` -o itpptest itpptest.cpp `itpp-config --libs`
    ```
    After this step, there should be an executable file,`itpptest` in itpptest folder
    ```shell
    $ ./itpptest
    ```
    You should see:
    ```
    Hello it++ library
    ```
- #### *Whole Ubuntu 16.04*
    If you installed Ubuntu on your old spare laptop, it is easier to install itpp library.
    - Install  **Synaptic Package Manager** onto your Ubuntu
    - Search for *\*it++ \**
    - Install all the packages from the results.

    Or, you can follow the instructions from this [video](https://www.youtube.com/watch?v=GWoVivaLzIo&t=199s)  

After installing ITPP on your Linux, we should compile the source file from `SWSC.zip`
In SWSC folder, there is one called *update*, under this folder:
```shell
$ vim makefile
```
Type in:
```makefile
SWSC:
    g++ -c main_swsc_single.cpp
    g++ -c turbo_rev.cpp
    g++ -c nit_coding.cpp
    g++ `itpp-config --cflags` -o program main_swsc_single.o turbo_rev.o nit_coding.o `itpp-config --libs`
```
Next, in command line:
```shell
$ make SWSC
$ ls
```
You should see a file named `program`, which is an executable file made by makefile.
Then,
```shell
$ ./program 512 1024 6 4.6 0.2 4000
```
After about 2 min, the result of this simulation will be shown.

#### Pipe line of C++ source files

Generate LDPC code ---> Treat Interference as Noise ---> P2P SWCM ---> 2-user SWSC vs. IAN

##### Specification of C++ files

- ldpc_gen.cpp

  This file is an example on [IT++ sourceforge](http://itpp.sourceforge.net/4.3.1/ldpc_gen_codes.html).

  It can generate different LDPC code and store them into .it files.

  The rate of these LDPC codes are all 0.5, but the component code rate may differ. Also, some of them are regular, others are irregular.

  ```shell
  $ make ldpc_gen
  $ ./ldpc_gen
  ```

- ldpc_gen_regular

  Generating customized LDPC *regular*code.

  **e.g.** To generate a LDPC code with code *Near = 1024*, component convolutional code *K = 3*, *N = 6*.

  ```shell
  $ make ldpc_gen_regular
  $ ./ldpc_gen_regular 1024 3 6
  ```

- IAN

  - point 2 point process showing below:

    message --> LDPC encoding --> modulation --> AWGN --> demodulation --> LDPC decoding

  - 2 user transmit messages following the above process, and one interferes the other.

  

  **e.g.** One code (2048_3_6) interfered by another (2048_3_8) and INR is 10dB.

  ```shell
  $ make IAN
  $ ./IAN LDPC_code_2048_3_6_R.it LDPC_code_2048_3_8_R.it 10
  ```

- SWCM_p2p

  message --> LDPC encoding --> SWCM encoding --> modulation --> AWGN --> demodulation --> LDPC & SWSC decoding

  **e.g.** Transfer one message encoded by LDPC(2048_3_6) and SWCM through AWGN, SNR = 10dB, block_num = 20

  ```sh
  $ make SWSC_p2p
  $ ./SWSC_p2p LDPC_code_2048_3_6_R.it 10 20
  ```

- SWSC_IAN

  One message transferred through SWCM channel, the other normal transmission would interfere SWCM channel.

  **e.g.** one LDPC code (2048_3_6) go through SWCM channel and another code(1024_3_6) interferes the former. (Ratio of Nvar MUST be 2:1)

  ```sh
  $ make SWSC_IAN
  $ ./SWSC_IAN LDPC_code_2048_3_6_R.it LDPC_code_1024_3_6_R.it 8 20
  ```

##### Make file

## Matlab source files

Most of the 

## UBC Sockeye Server