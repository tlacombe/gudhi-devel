CXX= g++-4.9 -O0 -gdwarf-3
CFLAGS= -std=c++11 -Wall -Wextra 

all: bottleneck.out

bottleneck.out: main.cpp
	$(CXX) $(CFLAGS) main.cpp -o bottleneck.out

clean:
	rm -f *.o bottleneck.out
