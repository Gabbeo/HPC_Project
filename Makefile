CC = cc
SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
BIN = ./bin
TARGET = heat_equation.out

LDFLAGS = -fopenmp
CFLAGS = -g -Wall -O3 -fopenmp -std=c11

all: dir $(BIN)/$(TARGET)

dir: ${BIN}

${BIN}:
	mkdir -p $(BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(BIN)/$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ -lm

.PHONY: clean
clean:
	rm -f $(OBJ) $(BIN)/$(TARGET)
