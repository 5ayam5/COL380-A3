CC=mpic++
CFLAGS=-std=c++17 -g -Wall -Wextra -O3 -fopenmp

run:
	mpirun --bind-to none ./main $(data_dir) $(top_k) $(user) $(out_file)

main: main.cpp
	$(CC) $(CFLAGS) main.cpp -o main

clean:
	rm -rf main $(out_dir) $(out_file)