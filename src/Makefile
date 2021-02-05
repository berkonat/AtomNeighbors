CC=gcc-10
CXX = g++-10
RM=rm -f
CPPFLAGS=
LDFLAGS=
LDLIBS=

FLAGS = -fopenmp -arch x86_64 -Wall -Werror -std=c++11 -lstdc++
LINK = -fopenmp -arch x86_64 -Wall -Werror -std=c++11 -lstdc++

all: AtomNeighs.x

OBJS= utils.o neigh.o main.o

AtomNeighs.x : $(OBJS)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) $(LINK) $(OBJS) -o $@

%.o : %.cpp
	$(CXX) $(CPPFLAGS) $(FLAGS) -c $<

clean :
	$(RM) *.o *.x
