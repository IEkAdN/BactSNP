#CXXFLAGS = -std=c++0x -MMD -O0 -Wall -g
CXXFLAGS = -std=c++0x -MMD -O3 -Wall

get_snp_phenix : main.o hoge.o fuga.o
	$(CXX) -o $@ $^

-include *.d
