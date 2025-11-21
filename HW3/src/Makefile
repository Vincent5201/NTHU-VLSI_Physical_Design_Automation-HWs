CC = g++
FLAGS = -Wall -O2

SRC = main.cpp deal_datas.cpp utils.cpp algos.cpp
TARGET = ../bin/hw3

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(FLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f ../bin/*