.PHONY: all clean

#TARGET = ./dphcar
TARGET1 = ./dph
TARGET2 = ./test
CC = gcc
CFLAGS = -Wall -Wextra -g -O0
LDFLAGS = -lm -lmpfr -lgmp
OBJS = dp2d.o fp.o globals.o rule.o

all: $(TARGET) $(TARGET1)

#$(TARGET): $(OBJS)

$(TARGET1): $(OBJS)

$(TARGET2): $(OBJS)

clean:
	@$(RM) $(OBJS) $(TARGET) $(TARGET1) $(TARGET2)
