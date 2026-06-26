CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libsofe.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_self_organization.c
TEST_BIN = build/test_self_organization
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))

.PHONY: all builddir clean test examples demo

all: builddir $(LIB)

builddir:
	@mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lsofe $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lsofe $(LDFLAGS)

demo: examples
	@echo "=== Brusselator Example ==="
	./build/example_brusselator
	@echo ""
	@echo "=== Turing Pattern Example ==="
	./build/example_turing
	@echo ""
	@echo "=== Benard Convection Example ==="
	./build/example_benard

clean:
	rm -rf build
