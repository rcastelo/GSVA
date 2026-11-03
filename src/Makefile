-include Makeconf

SHLIB = GSVA.so
OBJECTS = $(CU_OBJS) $(C_OBJS)

all: $(SHLIB)

$(SHLIB): $(OBJECTS)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(ALL_LIBS)

cuda/%.o: cuda/%.cu
	$(NVCC) -std=c++17 $(FP32_FLAG) $(GPU_ARCH_FLAGS) $(NVCC_FLAGS) -MMD -MP -c $< -o $@ $(HAVE_CUDA_DEF) $(ALL_CPPFLAGS) -I. -Xcompiler "$(CPICFLAGS) $(CFLAGS)"

-include cuda/*.d

clean:
	rm -f *.o cuda/*.o cuda/*.d $(SHLIB)

.PHONY: all clean
