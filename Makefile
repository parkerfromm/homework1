CPP      = g++  	# CPP compiler name
CPPFLAGS = -g 

EXCS = homework1

all: 
	make $(EXCS)
 
homework1: homework1.cpp 
	$(CPP) $(CPPFLAGS) homework1.cpp -o homework1

run:
	./$(EXCS)

clean:
	rm -f $(EXCS)