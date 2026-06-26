CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/liboep.a
SRCS = src/oep_core.c src/oep_order_parameter.c src/oep_slaving_principle.c \
       src/oep_phase_transition.c src/oep_emergence_metrics.c \
       src/oep_collective_dynamics.c src/oep_landau_ginzburg.c \
       src/oep_applications.c
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_oep.c
TEST_BIN = build/test_oep
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))

.PHONY: all clean test examples

all: builddir $(LIB)

builddir:
	mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -loep $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -loep $(LDFLAGS)

clean:
	rm -rf build
