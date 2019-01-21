#CXXFLAGS = -std=c++0x -MMD -O0 -Wall -g
CXXFLAGS = -std=c++0x -MMD -O3 -Wall

a : main.o hoge.o fuga.o
	$(CXX) -o $@ $^

-include *.d
