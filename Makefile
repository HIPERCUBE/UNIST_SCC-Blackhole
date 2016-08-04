CC = mpicc
CFLAGS = 
TARGET = blackhole_lab
OBJECTS = main.o blackhole_lab.o bundle.o KBLspacetime.o

all : $(TARGET)

$(TARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ -lm -lpng

clean :
	rm *.o $(TARGET)
