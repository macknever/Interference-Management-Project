ldpc_gen:
	g++ `itpp-config --cflags` -o ldpc_gen ldpc_gen.cpp `itpp-config --libs`

ldpc_simulation:
	g++ `itpp-config --cflags` -o ldpc_simulation LDPC_simulation.cpp `itpp-config --libs`

ldpc_gen_test:
	g++ `itpp-config --cflags` -o ldpc_gen_test ldpc_gen_test.cpp `itpp-config --libs`

ldpc_simu_test:
	rm -rf ldpc_simu_test
	g++ `itpp-config --cflags` -o ldpc_simu_test ldpc_simu_test.cpp `itpp-config --libs`

ldpc_gen_regular:
	g++ `itpp-config --cflags` -o ldpc_gen_regular ldpc_gen_regular.cpp `itpp-config --libs`
test:
	
	g++ `itpp-config --cflags` -o test test.cpp `itpp-config --libs`

SWSC_IAN:
	g++ -c -g SWSC.cpp
	g++ -c -g SWSC_LDPC_SIMU.cpp
	g++ `itpp-config --cflags` -o SWSC_IAN SWSC.o SWSC_LDPC_SIMU.o `itpp-config --libs`
SWSC_p2p:
	g++ -c -g SWSC.cpp
	g++ -c -g SWSC_p2p.cpp
	g++ `itpp-config --cflags` -o SWSC_p2p SWSC.o SWSC_p2p.o `itpp-config --libs`
IAN:
	g++ `itpp-config --cflags` -o IAN IAN_LDPC_SIMU.cpp `itpp-config --libs`

# SWSC:
# 	g++ `itpp-config --cflags` -o SWSC SWSC.cpp `itpp-config --libs`