CXX      = g++
CXXFLAGS = -s -std=c++0x -O3 -funroll-loops -fomit-frame-pointer -fopenmp
# CXXFLAGS = -s -std=c++0x -O0 -funroll-loops -fopenmp -g3
LDFLAGS  = 

INTERNAL_DEF = -DRUN_MODE=\"internal\"
TRIM_DEF = -DRUN_MODE=\"\"
#OBJS = common.o dense.o ecc.o
OBJS_TRIM = common.o trim.o main.o
OBJS_INTERNAL = common.o trim.o main.oo
TARGET = platanus_trim platanus_internal_trim


all: $(TARGET)

platanus_trim: $(OBJS_TRIM)
	$(CXX) $(CXXFLAGS) -o $@ $^

platanus_internal_trim: $(OBJS_INTERNAL)
	$(CXX) $(CXXFLAGS) -o $@ $^
.cpp.o:
	$(CXX) -o $@ -c $< $(CXXFLAGS) $(TRIM_DEF)
main.oo: main.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS) $(INTERNAL_DEF)

clean:
	rm -f $(OBJS_TRIM) $(TARGET) $(OBJS_INTERNAL)

