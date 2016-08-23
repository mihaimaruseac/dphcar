.PHONY: all clean

TARGET = ./dph
CC = gcc
CFLAGS = -Wall -Wextra -g -O0
LDFLAGS = -lm
OBJS = dp2d.o fp.o globals.o rule.o histogram.o rs.o

all: $(TARGET)

$(TARGET): $(OBJS)

clean:
	@$(RM) $(OBJS) $(TARGET)
