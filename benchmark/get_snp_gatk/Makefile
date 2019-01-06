#CXXFLAGS = -std=c++0x -MMD -O0 -Wall -g
CXXFLAGS = -std=c++0x -MMD -O3 -Wall

get_snp_gatk : main.o two_matrix.o
	$(CXX) -o $@ $^

-include *.d
