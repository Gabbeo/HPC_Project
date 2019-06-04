CC = mpicc
SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
BIN = ./bin
TARGET = heat_equation.out

LDFLAGS = -lm  -fopenmp
CFLAGS = -g -Wall -O3 -fopenmp

all: dir $(BIN)/$(TARGET)

dir: ${BIN}

${BIN}:
	mkdir -p $(BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(BIN)/$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm -f $(OBJ) $(BIN)/$(TARGET)
