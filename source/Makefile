CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS = Conn.o CSize.o

LIBS =

TARGET = Conn CSize

Conn: Conn.o
	$(CXX) -o Conn Conn.o

CSize: CSize.o
	$(CXX) -o CSize CSize.o

all: $(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)


