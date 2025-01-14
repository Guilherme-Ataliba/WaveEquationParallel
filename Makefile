# Compilers
CC = gcc
CUDA_CC = nvcc
MPI_CC = mpicc
PY = python

# Flags
CFLAGS = -lm 
CUDA_FLAGS = 
MPI_FLAGS = -lm
OMP_FLAGS = -fopenmp -lm

# Directories
INCLUDE = include
BUILD = build
SRC = source
BIN = scripts
FIG = figures
OUTPUT = output
MULTI_MATRIX = $(OUTPUT)/multi_matrix  # Directory used to house each matrix output

# Files
GIF_FILE = figure.gif

# Dependencies
DEPS = $(INCLUDE)/auxilary.h
SERIAL_OBJS = $(BUILD)/algorithm_serial.o $(BUILD)/auxilary_serial.o
CUDA_OBJS = $(BUILD)/algorithm_cuda.o $(BUILD)/auxilary_cuda.o
MPI_OBJS = $(BUILD)/algorithm_mpi.o $(BUILD)/auxilary_mpi.o
OMP_OBJS = $(BUILD)/algorithm_omp.o $(BUILD)/auxilary_omp.o

# Target to choose between serial, CUDA, MPI, and OpenMP
USE_CUDA = 0  # Set to 1 for CUDA, 0 for Serial
MPI_CORES ?= 4  # Default number of cores for MPI
OMP_THREADS ?= 4  # Default number of threads for OpenMP

# Compiler targets
.PHONY: all
all: $(SRC)/algorithm_serial

# Run target for Serial
.PHONY: run_serial
run_serial: $(SRC)/algorithm_serial | $(SRC) $(BIN) $(OUTPUT) $(MULTI_MATRIX)
	./$(BIN)/algorithm_serial
	make clean

# Run target for CUDA
.PHONY: run_cuda
run_cuda: $(SRC)/algorithm_cuda | $(SRC) $(BIN) $(OUTPUT) $(MULTI_MATRIX)
	./$(BIN)/algorithm_cuda
	make clean

# Run target for MPI
.PHONY: run_mpi
run_mpi: $(SRC)/algorithm_mpi | $(SRC) $(BIN) $(OUTPUT) $(MULTI_MATRIX)
	mpirun -np $(MPI_CORES) ./$(BIN)/algorithm_mpi
	make clean

# Run target for OpenMP
.PHONY: run_omp
run_omp: $(SRC)/algorithm_omp | $(SRC) $(BIN) $(OUTPUT) $(MULTI_MATRIX)
	OMP_NUM_THREADS=$(OMP_THREADS) ./$(BIN)/algorithm_omp
	make clean

# Plot target
.PHONY: plot
plot: 
	$(PY) $(SRC)/plot.py

# Compile Serial Algorithm
$(SRC)/algorithm_serial: $(SERIAL_OBJS) | $(SRC) $(BIN)
	$(CC) $(SERIAL_OBJS) -o $(BIN)/$(@F) $(CFLAGS)

# Compile CUDA Algorithm
$(SRC)/algorithm_cuda: $(CUDA_OBJS) | $(SRC) $(BIN)
	$(CUDA_CC) $(CUDA_OBJS) -o $(BIN)/$(@F) $(CUDA_FLAGS)

# Compile MPI Algorithm
$(SRC)/algorithm_mpi: $(MPI_OBJS) | $(SRC) $(BIN)
	$(MPI_CC) $(MPI_OBJS) -o $(BIN)/$(@F) $(MPI_FLAGS)

# Compile OpenMP Algorithm
$(SRC)/algorithm_omp: $(OMP_OBJS) | $(SRC) $(BIN)
	$(CC) $(OMP_OBJS) -o $(BIN)/$(@F) $(OMP_FLAGS)

# Object files for Serial Code
$(BUILD)/algorithm_serial.o: $(SRC)/algorithm_serial.c $(DEPS) | $(SRC) $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/auxilary_serial.o: $(SRC)/auxilary.c $(DEPS) | $(SRC) $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

# Object files for CUDA Code
$(BUILD)/algorithm_cuda.o: $(SRC)/algorithm_cuda.cu $(DEPS) | $(SRC) $(BUILD)
	$(CUDA_CC) $(CUDA_FLAGS) -c $< -o $@

$(BUILD)/auxilary_cuda.o: $(SRC)/auxilary.c $(DEPS) | $(SRC) $(BUILD)
	$(CUDA_CC) $(CUDA_FLAGS) -c $< -o $@

# Object files for MPI Code
$(BUILD)/algorithm_mpi.o: $(SRC)/algorithm_mpi.c $(DEPS) | $(SRC) $(BUILD)
	$(MPI_CC) $(MPI_FLAGS) -c $< -o $@

$(BUILD)/auxilary_mpi.o: $(SRC)/auxilary.c $(DEPS) | $(SRC) $(BUILD)
	$(MPI_CC) $(MPI_FLAGS) -c $< -o $@

# Object files for OpenMP Code
$(BUILD)/algorithm_omp.o: $(SRC)/algorithm_omp.c $(DEPS) | $(SRC) $(BUILD)
	$(CC) $(OMP_FLAGS) -c $< -o $@

$(BUILD)/auxilary_omp.o: $(SRC)/auxilary.c $(DEPS) | $(SRC) $(BUILD)
	$(CC) $(OMP_FLAGS) -c $< -o $@

# Clean target
.PHONY: clean
clean: | $(BUILD)
	rm -f $(BUILD)/*.o

# Remove target
.PHONY: remove
remove: | $(BIN)
	rm -f $(BIN)/*.exe
	rm -f $(OUTPUT)/Execution_time.txt
	rm -f $(OUTPUT)/multi_matrix/*

# Graph target (GIF)
.PHONY: graph
graph: $(FIG)/$(GIF_FILE) | $(FIG)
	start "$(FIG)/$(GIF_FILE)"

# Make the GIF
$(FIG)/$(GIF_FILE): | $(FIG)
	make run_serial

# Directories Creation
$(INCLUDE):
	mkdir $@

$(BUILD):
	mkdir $@

$(SRC):
	mkdir $@

$(BIN):
	mkdir $@

$(FIG):
	mkdir $@

$(OUTPUT):
	mkdir $@

$(MULTI_MATRIX): 
	mkdir $@
