CC       = gcc
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra -Wpedantic -Wshadow
LDLIBS  ?= -lm

ifeq ($(OS),Windows_NT)
EXE = .exe
else
EXE =
endif

BIN      = project1$(EXE)
SRC      = src/main.c src/cluster.c src/diffusion.c
OBJ      = $(SRC:.c=.o)

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDLIBS)

src/%.o: src/%.c src/cluster.h src/diffusion.h
	$(CC) $(CFLAGS) -c -o $@ $<

check: $(BIN)
	./$(BIN) --check

conv: $(BIN)
	./$(BIN) --conv

figures:
	cd tools && node plot.mjs

clean:
	$(RM) $(OBJ) $(BIN) project1 project1.exe
	$(RM) -r output

.PHONY: all check conv figures clean
