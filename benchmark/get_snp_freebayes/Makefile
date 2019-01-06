#CXXFLAGS = -std=c++0x -MMD -O0 -Wall -g
CXXFLAGS = -std=c++0x -MMD -O3 -Wall

get_snp_freebayes : main.o hoge.o
	$(CXX) -o $@ $^

-include *.d
