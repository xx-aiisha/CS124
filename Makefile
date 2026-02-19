CXX = clang++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -pedantic

randmst: randmst.cpp
	$(CXX) $(CXXFLAGS) randmst.cpp -o randmst

