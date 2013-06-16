# Makefile

CC    =icc
DEVCC =nvcc
LINK =nvcc 
LIB        =-ccbin=icc -Xcompiler "-mkl:parallel -openmp -O3"
DEVLIB = -lcublas -lcusparse
CFLAGS          =-c -O3 
CFLAGSDEBUG     =-c -g -G -O0
SRC_DIR         = src
INCLUDE_DIR     =-I $(SRC_DIR)
BUILD_DIR       = build
RESULT_DIR      = .
NAME            = solver_mkl_cuda.exe
BASE_INCLUDE_DIR= $(BASE_SRC_DIR)
OPTIMAZE_SPECIFIC =-openmp -msse4.1
# Linker flags
LDFLAGS         = $(DEVLIB) $(LIB) $(INCLUDE_DIR)



RESULT = $(RESULT_DIR)/$(NAME)

OBJS = \
        $(BUILD_DIR)/main.o                     \
        $(BUILD_DIR)/bilu0.o               \
        $(BUILD_DIR)/mmio.o                 \
        $(BUILD_DIR)/solver_cusparse.o               \
	$(BUILD_DIR)/solver_mkl.o               \
	$(BUILD_DIR)/utils.o               \
	$(BUILD_DIR)/dev_utils.o
	
################################################################################
# Build rules
all : $(RESULT)

$(RESULT) : $(OBJS)
	$(LINK) $(LDFLAGS) $(OBJS) -o $(RESULT)

clean :
	rm -f $(OBJS) $(RESULT)

############################################################################
# Individual files build rules
$(BUILD_DIR)/main.o        : $(SRC_DIR)/main.cpp
	$(DEVCC) $(CFLAGS) -ccbin=$(CC) -Xcompiler "$(OPTIMAZE_SPECIFIC)" $< -o $@
$(BUILD_DIR)/dev_utils.o        : $(SRC_DIR)/dev_utils.cu
	$(DEVCC) $(CFLAGS) $< -o $@
$(BUILD_DIR)/solver_cusparse.o        : $(SRC_DIR)/solver_cusparse.cu
	$(DEVCC) $(CFLAGS) -ccbin=icc -Xcompiler "-openmp" $< -o $@
#	$(DEVCC) $(CFLAGS) -ccbin=icc $< -o $@
$(BUILD_DIR)/bilu0.o        : $(SRC_DIR)/bilu0.cpp
	$(CC) $(CFLAGS) $(OPTIMAZE_SPECIFIC) $< -o $@
$(BUILD_DIR)/mmio.o        : $(SRC_DIR)/mmio.cpp
	$(CC) $(CFLAGS) $< -o $@
$(BUILD_DIR)/utils.o        : $(SRC_DIR)/utils.cpp
	$(CC) $(CFLAGS) $< -o $@
$(BUILD_DIR)/solver_mkl.o        : $(SRC_DIR)/solver_mkl.cpp
	$(CC) $(CFLAGS) $(OPTIMAZE_SPECIFIC) $< -o $@
