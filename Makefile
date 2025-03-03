# Compiler and flags
FC = gfortran
FFLAGS_BASE = -cpp -fopenmp

# Maximum optimization for cluster builds
FFLAGS_CLUSTER = $(FFLAGS_BASE) \
                 -Ofast \
                 -march=native \
                 -mtune=native \
                 -ffast-math \
                 -funroll-loops \
                 -ftree-vectorize \
                 -fno-stack-arrays \
                 -flto \
                 -fno-range-check \
                 -fno-signed-zeros \
                 -fno-trapping-math \
                 -mfpmath=sse \
                 -malign-data=cacheline \
                 -fipa-pta \
                 -floop-nest-optimize

# Local release optimization (still fast but more portable)
FFLAGS_OPT = $(FFLAGS_BASE) \
             -O3 \
             -march=native \
             -ffast-math \
             -funroll-loops \
             -ftree-vectorize \
             -flto

# Debug flags
FFLAGS_DEBUG = $(FFLAGS_BASE) \
               -g \
               -fcheck=all \
               -fbacktrace \
               -Wall \
               -Wextra \
               -Warray-temporaries \
               -fimplicit-none \
               -finit-real=nan \
               -finit-integer=-99999999

# Profile flags
FFLAGS_PROF = $(FFLAGS_BASE) -pg -g

# Perf flags
FFLAGS_PERF = $(FFLAGS_BASE) -O2 -g

# Source files in dependency order
SRCS = types_mod.f90 \
       error_mod.f90 \
       config_mod.f90 \
       cell_list_mod.f90 \
       coordination_mod.f90 \
       progress_mod.f90 \
       data_file_io_mod.f90 \
       trajectory_io_mod.f90 \
       output_io_mod.f90 \
       io_mod.f90 \
       benchmark_mod.f90 \
       main.f90

# Object files
OBJS = $(SRCS:.f90=.o)

# Module files
MODS = $(SRCS:.f90=.mod)

# Executables
PROG = coordination_analysis
PROG_CLUSTER = $(PROG)_cluster
PROG_DEBUG = $(PROG)_debug
PROG_PROF = $(PROG)_prof
PROG_PERF = $(PROG)_perf

# Default target
all: release cluster debug profile perf

# Cluster build (maximum optimization)
cluster: FFLAGS = $(FFLAGS_CLUSTER)
cluster: $(PROG_CLUSTER)

# Release build (local optimized)
release: FFLAGS = $(FFLAGS_OPT)
release: $(PROG)

# Debug build
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: $(PROG_DEBUG)

# Profile build
profile: FFLAGS = $(FFLAGS_PROF)
profile: $(PROG_PROF)

# Perf build
perf: FFLAGS = $(FFLAGS_PERF)
perf: $(PROG_PERF)

# Build rules
$(PROG): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(PROG_CLUSTER): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(PROG_DEBUG): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(PROG_PROF): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(PROG_PERF): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

# Call graph target
callgraph: $(PROG_PROF)
	./$(PROG_PROF)
	gprof -q -b $(PROG_PROF) gmon.out > callgraph.txt

# Dependencies
types_mod.o: types_mod.f90
error_mod.o: error_mod.f90
config_mod.o: config_mod.f90 types_mod.o error_mod.o
cell_list_mod.o: cell_list_mod.f90 types_mod.o config_mod.o error_mod.o
coordination_mod.o: coordination_mod.f90 types_mod.o config_mod.o cell_list_mod.o error_mod.o
progress_mod.o: progress_mod.f90 types_mod.o config_mod.o
data_file_io_mod.o: data_file_io_mod.f90 types_mod.o config_mod.o error_mod.o
trajectory_io_mod.o: trajectory_io_mod.f90 types_mod.o config_mod.o error_mod.o
output_io_mod.o: output_io_mod.f90 types_mod.o config_mod.o error_mod.o coordination_mod.o cell_list_mod.o
io_mod.o: io_mod.f90 types_mod.o trajectory_io_mod.o data_file_io_mod.o output_io_mod.o progress_mod.o
benchmark_mod.o: benchmark_mod.f90 types_mod.o config_mod.o cell_list_mod.o coordination_mod.o io_mod.o error_mod.o
main.o: main.f90 types_mod.o config_mod.o cell_list_mod.o coordination_mod.o io_mod.o error_mod.o benchmark_mod.o

# Compilation rule
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Clean rules
clean:
	rm -f *.o *.mod $(PROG) $(PROG_CLUSTER) $(PROG_DEBUG) $(PROG_PROF) $(PROG_PERF)

clean_prof:
	rm -f gmon.out profile.txt callgraph.txt

clean_all: clean clean_prof

# Help target
help:
	@echo "Available targets:"
	@echo "  make        - Build all versions"
	@echo "  make cluster - Build highly optimized version for cluster"
	@echo "  make callgraph - Generate call graph using gprof"
	@echo "  make release - Build optimized version for local use"
	@echo "  make debug   - Build debug version"
	@echo "  make profile - Build version with profiling"
	@echo "  make perf    - Build version with perf profiling"
	@echo "  make clean   - Remove compiled files"
	@echo "  make clean_prof - Remove profiling output"
	@echo "  make clean_all  - Remove all generated files"
	@echo "  make help    - Show this help message"

.PHONY: all cluster release debug profile perf clean clean_prof clean_all help