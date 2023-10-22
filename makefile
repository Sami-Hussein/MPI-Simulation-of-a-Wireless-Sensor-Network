CC = mpicc
CFLAGS = -lm -fopenmp -Wall
TARGET = test

all: $(TARGET)

$(TARGET): wsn.c
	$(CC) wsn.c -o $(TARGET) $(CFLAGS)

run1: $(TARGET)
	mpirun -oversubscribe -np 17 ./$(TARGET) 4 4 20

run2: $(TARGET)
	mpirun -oversubscribe -np 10 ./$(TARGET) 3 3 20

clean:
	rm -f $(TARGET)
