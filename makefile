all:
	g++ -std=c++11  -O3 -o solver src/solver.cpp -fopenmp

clean:
	$(RM) solver
