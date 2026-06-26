CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libprigogine.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRCS = $(wildcard tests/*.c)
TEST_BINS = $(patsubst tests/%.c,build/%,$(TEST_SRCS))
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))

.PHONY: all builddir clean test examples demo

all: builddir $(LIB)

builddir:
	mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: all $(TEST_BINS)
	for t in $(TEST_BINS); do ./$$t; done

build/%: tests/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lprigogine $(LDFLAGS)

examples: all $(EXAMPLE_BINS)

build/%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lprigogine $(LDFLAGS)

demo: examples
	for e in $(EXAMPLE_BINS); do echo "=== $$e ==="; ./$$e; done

clean:
	rm -rf build
