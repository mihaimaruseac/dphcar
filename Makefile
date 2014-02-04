.PHONY: all clean

TARGET = ./dphcar
CC = gcc
CFLAGS = -Wall -Wextra -g -O0
LDFLAGS = -lm -lmpfr -lgmp
OBJS = 

all: $(TARGET)

$(TARGET): $(OBJS)

clean:
	@$(RM) $(OBJS) $(TARGET)
