# Makefile for mini-entropy-production-minimization
# Builds the entropy production minimization library and tests
#
# Usage:
#   make          - Build library and tests
#   make test     - Build and run tests
#   make examples - Build example programs
#   make clean    - Remove build artifacts

CC       = gcc
CFLAGS   = -Wall -Wextra -std=c99 -O2 -g
LDFLAGS  = -lm

# Directories
INCDIR   = include
SRCDIR   = src
TESTDIR  = tests
EXAMPLEDIR = examples
DOCSDIR  = docs
BUILDDIR = build

# Source files
SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(SRCS:$(SRCDIR)/%.c=$(BUILDDIR)/%.o)

# Header files
HDRS = $(wildcard $(INCDIR)/*.h)

# Test and example targets
TESTS = $(wildcard $(TESTDIR)/*.c)
TESTBINS = $(TESTS:$(TESTDIR)/%.c=$(BUILDDIR)/test_%)

EXAMPLES = $(wildcard $(EXAMPLEDIR)/*.c)
EXAMPLEBINS = $(EXAMPLES:$(EXAMPLEDIR)/%.c=$(BUILDDIR)/example_%)

.PHONY: all clean test examples

all: $(BUILDDIR)/libepm.a

# Build directory
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Compile sources to objects
$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(HDRS) | $(BUILDDIR)
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $@

# Build static library
$(BUILDDIR)/libepm.a: $(OBJS)
	ar rcs $@ $(OBJS)

# Build and run tests
test: $(TESTBINS)
	@echo "=== Running tests ==="
	@passed=0; failed=0; \
	for t in $(TESTBINS); do \
		echo "Running $$t..."; \
		if ./$$t; then \
			passed=$$((passed+1)); \
		else \
			failed=$$((failed+1)); \
		fi; \
	done; \
	echo "=== Results: $$passed passed, $$failed failed ==="; \
	if [ $$failed -gt 0 ]; then exit 1; fi

# Link test binaries
$(BUILDDIR)/test_%: $(TESTDIR)/%.c $(BUILDDIR)/libepm.a | $(BUILDDIR)
	$(CC) $(CFLAGS) -I$(INCDIR) $< -L$(BUILDDIR) -lepm $(LDFLAGS) -o $@

# Build examples
examples: $(EXAMPLEBINS)

$(BUILDDIR)/example_%: $(EXAMPLEDIR)/%.c $(BUILDDIR)/libepm.a | $(BUILDDIR)
	$(CC) $(CFLAGS) -I$(INCDIR) $< -L$(BUILDDIR) -lepm $(LDFLAGS) -o $@

# Clean
clean:
	rm -rf $(BUILDDIR)

# Line count
count:
	@echo "Header files:"
	@wc -l $(INCDIR)/*.h
	@echo ""
	@echo "Source files:"
	@wc -l $(SRCDIR)/*.c
	@echo ""
	@echo -n "Total: "
	@cat $(INCDIR)/*.h $(SRCDIR)/*.c | wc -l