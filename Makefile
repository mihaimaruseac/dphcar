.PHONY: all clean

TARGET = ./dph
NPTARGET = ./nph
CC = gcc
CFLAGS = -Wall -Wextra -O3 #g -O0
LDFLAGS = -lm -lmpfr -lgmp
OBJS = dp2d.o fp.o globals.o rule.o histogram.o

all: $(TARGET) $(NPTARGET)

$(TARGET): $(OBJS)

$(NPTARGET): $(OBJS)

clean:
	@$(RM) $(OBJS) $(TARGET) $(NPTARGET)
