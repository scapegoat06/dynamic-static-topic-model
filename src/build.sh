g++ -Wall -O3 -mtune=native -march=native -c dstm.cc
g++ -Wall -O3 -mtune=native -march=native -c main.cc
g++ -Wall -O3 -mtune=native -march=native -o dstm dstm.o dstm.o
rm dstm.o dstm.o
