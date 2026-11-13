-include Makeconf

SHLIB = GSVA.so
OBJECTS = $(CU_OBJS) $(C_OBJS)

all: $(SHLIB)

$(SHLIB): $(OBJECTS)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(ALL_LIBS)

cuda_params.o: cuda_params.c cuda/config.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(FP32_FLAG) -c $< -o $@

cuda/%.o: cuda/%.cu
	$(NVCC) -std=c++17 $(GPU_ARCH_FLAGS) $(NVCC_FLAGS) -MMD -MP -c $< -o $@ $(ALL_CPPFLAGS) -I. -Xcompiler "$(CPICFLAGS) $(CFLAGS)"

-include cuda/*.d

clean:
	rm -f *.o cuda/*.o cuda/*.d $(SHLIB)

.PHONY: all clean
