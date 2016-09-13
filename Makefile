.PHONY: all clean

TARGET = ./dph ./cr
CC = gcc
CFLAGS = -Wall -Wextra -g -O0
LDFLAGS = -lm
OBJS = rs.o fp.o globals.o histogram.o itstree.o recall.o dp2d.o

all: $(TARGET)

$(TARGET): $(OBJS)

clean:
	@$(RM) $(OBJS) $(TARGET)
