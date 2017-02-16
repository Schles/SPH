X = g++
INCLUDES= -I./../../extern/eigen/include
CXXFLAGS = -g $(INCLUDES)
SRCM= ../Kernel.cpp
OBJM = $(SRCM:.cpp=.o)
	LINKFLAGS= -lcppunit

kerneltest: KernelTest.cpp $(OBJM)
		$(CXX) $(CXXFLAGS) -o $@ KernelTest.cpp $(OBJM) $(LINKFLAGS) $(LINKFLAGSLOG4) $(LIBLOG)

# Default compile
#
# .cpp.o:
# 	$(CXX) $(CXXFLAGS) -c $< -o $@
