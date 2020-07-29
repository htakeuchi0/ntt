# compilers
CXX = g++
# CXXFLAGS = -std=c++17 -O3 -g -Wall --pedantic-errors -fsanitize=address -fno-omit-frame-pointer
CXXFLAGS = -std=c++17 -O3 -g -Wall
MAIN_TARGET = build
TEST_TARGET = test

# sources
SRCS = src/*.cpp
MAIN_SRCS = main/main.cpp $(SRCS)
TEST_SRCS = test/*.cpp $(SRCS)

# includes
INCLUDES = -I. 
MAIN_INCLUDES = $(INCLUDES)
TEST_INCLUDES = $(INCLUDES) -Igoogletest-release-1.10.0/googletest/include

# for test
TEST_LIBS = -Lgoogletest-release-1.10.0/googletest/build/lib
TEST_LINKS = -lgtest -lgtest_main -lpthread

# objects
MAIN_OBJ = main.o
TEST_OBJ = test.o

# doxygen
DOXYGEN = doxygen
BROWSER = firefox
INDEXPATH = doxygen/html/index.html

.PHONY: all clean main test docs

all: $(MAIN_TARGET) $(TEST_TARGET)

$(MAIN_TARGET): $(MAIN_SRCS)
	$(CXX) $(CXXFLAGS) $(MAIN_SRCS) $(MAIN_INCLUDES) -o $(MAIN_OBJ)

$(TEST_TARGET): $(TESTS_SRCS)
	$(CXX) $(CXXFLAGS) $(TEST_SRCS) $(TEST_INCLUDES) $(TEST_LIBS) $(TEST_LINKS) -o $(TEST_OBJ)

clean:
	rm -f *.o

docs:
	$(DOXYGEN)
	$(BROWSER) $(INDEXPATH) &
