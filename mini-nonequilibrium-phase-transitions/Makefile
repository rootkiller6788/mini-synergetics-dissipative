CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libnept.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_nept.c
TEST_BIN = build/test_nept
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
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lnept $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lnept $(LDFLAGS)

demo: examples
	@echo "=== Running demos ==="
	@for exe in $(EXAMPLE_BINS); do echo "--- $$exe ---"; ./$$exe || true; done

clean:
	rm -rf build