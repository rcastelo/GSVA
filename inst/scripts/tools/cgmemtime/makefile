
CFLAGSW_GCC = -Wall -Wextra -Wno-missing-field-initializers \
    -Wno-missing-braces \
    -Wmissing-prototypes -Wfloat-equal \
    -Wwrite-strings -Wpointer-arith -Wcast-align \
    -Wnull-dereference \
    -Werror=multichar -Werror=sizeof-pointer-memaccess -Werror=return-type \
    -fstrict-aliasing

CFLAGS = -g -Og $(CFLAGSW_GCC)


EXEC = cgmemtime testa
TEMP = $(EXEC)


.PHONY: all
	
all: $(EXEC)

.PHONY: clean

clean:
	rm -f $(TEMP)


.PHONY: check

check:
	bash test.sh
