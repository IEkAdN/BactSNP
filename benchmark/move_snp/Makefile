CXXFLAGS = -std=c++0x -MMD -O0 -Wall -g
#CXXFLAGS = -std=c++0x -MMD -O3 -Wall

move_snp : main.o hoge.o fuga.o
	$(CXX) -o $@ $^

-include *.d
