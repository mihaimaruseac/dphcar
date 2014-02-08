.PHONY: all clean

TARGET = ./dphcar
TARGET1 = ./dph
CC = gcc
CFLAGS = -Wall -Wextra -g -O0
LDFLAGS = -lm -lmpfr -lgmp
OBJS = dp2d.o fp.o globals.o

all: $(TARGET) $(TARGET1)

$(TARGET): $(OBJS)

$(TARGET1): $(OBJS)

clean:
	@$(RM) $(OBJS) $(TARGET) $(TARGET1)
