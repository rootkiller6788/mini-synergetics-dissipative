CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libhaken.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRCS = $(wildcard tests/*.c)
TEST_BINS = $(patsubst tests/%.c,build/test_%,$(TEST_SRCS))
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/example_%,$(EXAMPLE_SRCS))
DEMO_SRCS = $(wildcard demos/*.c)
DEMO_BINS = $(patsubst demos/%.c,build/demo_%,$(DEMO_SRCS))
BENCH_SRCS = $(wildcard benches/*.c)
BENCH_BINS = $(patsubst benches/%.c,build/bench_%,$(BENCH_SRCS))
all: builddir $(LIB)
builddir: ; @mkdir -p build
$(LIB): $(OBJS) ; $(AR) rcs $@ $^
build/%.o: src/%.c | builddir ; $(CC) $(CFLAGS) -c $< -o $@
build/test_%: tests/%.c $(LIB) ; $(CC) $(CFLAGS) $< -Lbuild -lhaken $(LDFLAGS) -o $@
build/example_%: examples/%.c $(LIB) ; $(CC) $(CFLAGS) $< -Lbuild -lhaken $(LDFLAGS) -o $@
build/demo_%: demos/%.c $(LIB) ; $(CC) $(CFLAGS) $< -Lbuild -lhaken $(LDFLAGS) -o $@
build/bench_%: benches/%.c $(LIB) ; $(CC) $(CFLAGS) $< -Lbuild -lhaken $(LDFLAGS) -o $@
test: $(LIB) $(TEST_BINS) ; for t in $(TEST_BINS); do $$t || exit 1; done && echo All tests passed
examples: $(LIB) $(EXAMPLE_BINS) ; for e in $(EXAMPLE_BINS); do $$e; done
demos: $(LIB) $(DEMO_BINS) ; for d in $(DEMO_BINS); do $$d; done
bench: $(LIB) $(BENCH_BINS) ; for b in $(BENCH_BINS); do $$b; done
clean: ; rm -rf build
.PHONY: all test examples demos bench clean builddir
