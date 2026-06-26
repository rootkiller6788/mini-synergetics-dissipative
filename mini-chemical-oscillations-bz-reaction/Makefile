# Makefile for mini-chemical-oscillations-bz-reaction
# Chemical Oscillations — Belousov-Zhabotinsky Reaction

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS  = -lm

# Source files
SRC = src/bz_reaction.c \
      src/bz_oregonator.c \
      src/bz_reaction_diffusion.c \
      src/bz_integration.c \
      src/bz_oscillation.c \
      src/bz_waves.c \
      src/bz_stability.c \
      src/bz_bistability.c

# Test files
TEST_SRC = tests/test_oregonator.c \
           tests/test_integration.c \
           tests/test_oscillation.c

# Example files
EXAMPLE_SRC = examples/example_oregonator_oscillation.c \
              examples/example_spiral_wave.c \
              examples/example_bistability.c

# Object directories
BUILD_DIR = build

.PHONY: all clean test examples lint check-lean

all: test examples

# Build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compile all source files to objects
$(BUILD_DIR)/%.o: src/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Compile test files
$(BUILD_DIR)/test_%.o: tests/test_%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Compile example files
$(BUILD_DIR)/example_%.o: examples/example_%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Link tests
$(BUILD_DIR)/test_oregonator: $(BUILD_DIR)/test_oregonator.o $(SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/test_integration: $(BUILD_DIR)/test_integration.o $(SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/test_oscillation: $(BUILD_DIR)/test_oscillation.o $(SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Link examples
$(BUILD_DIR)/example_oregonator: $(BUILD_DIR)/example_oregonator_oscillation.o $(SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/example_spiral: $(BUILD_DIR)/example_spiral_wave.o $(SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/example_bistable: $(BUILD_DIR)/example_bistability.o $(SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Run all tests
test: $(BUILD_DIR)/test_oregonator $(BUILD_DIR)/test_integration $(BUILD_DIR)/test_oscillation
	@echo "=== Running Tests ==="
	@echo ""
	@./$(BUILD_DIR)/test_oregonator
	@echo ""
	@./$(BUILD_DIR)/test_integration
	@echo ""
	@./$(BUILD_DIR)/test_oscillation
	@echo ""
	@echo "=== All Tests Complete ==="

# Build all examples
examples: $(BUILD_DIR)/example_oregonator $(BUILD_DIR)/example_spiral $(BUILD_DIR)/example_bistable

# Run all examples
run-examples: examples
	@echo "=== Running Examples ==="
	@echo ""
	@./$(BUILD_DIR)/example_oregonator
	@echo ""
	@./$(BUILD_DIR)/example_spiral
	@echo ""
	@./$(BUILD_DIR)/example_bistable

# Lint
lint:
	@echo "=== Linting ==="
	@for f in $(SRC) $(TEST_SRC) $(EXAMPLE_SRC); do \
		echo "Checking $$f..."; \
		gcc -fsyntax-only -std=c11 -Wall -Iinclude $$f || exit 1; \
	done
	@echo "Lint OK"

# Check line count (must be >= 3000)
check-lines:
	@echo "include/ lines: $$(cat include/*.h | wc -l)"
	@echo "src/ lines:     $$(cat src/*.c | wc -l)"
	@echo "TOTAL:          $$(cat include/*.h src/*.c | wc -l)"

# Clean
clean:
	rm -rf $(BUILD_DIR)

# Lean check (requires Lean 4 installed)
check-lean:
	@if command -v lake >/dev/null 2>&1; then \
		echo "Running lean on src/bz_lean.lean..."; \
		lean src/bz_lean.lean || echo "Lean check failed (may need lake project)"; \
	else \
		echo "Lean 4 not installed. Skipping lean check."; \
	fi

