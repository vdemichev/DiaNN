.PHONY: all
all: test
# This rule tells make how to build hello from test1.cpp
test: test2.cpp
	g++ -O2 -I../../include test2.cpp -o test2.o

# This rule tells make to delete the program
.PHONY: clean 
clean:
	rm -f test1.o


