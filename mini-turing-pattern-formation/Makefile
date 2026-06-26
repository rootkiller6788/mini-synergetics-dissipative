CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libturing_pattern.a

SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))

TEST_SRC = tests/test_main.c
TEST_BIN = build/test_turing

EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))

.PHONY: all builddir clean test examples demo

all: builddir $(LIB)

builddir:
	@mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c builddir
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lturing_pattern $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/%: examples/%.c $(LIB) builddir
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lturing_pattern $(LDFLAGS)

build/example_%: examples/example_%.c $(LIB) builddir
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lturing_pattern $(LDFLAGS)

demo: examples
	@echo "=== Running example_stability_analysis ==="
	./build/example_stability_analysis

clean:
	rm -rf build
