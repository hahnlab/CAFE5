# compilers and flags
CC = g++

# DEBUG=true

ifdef DEBUG
  CFLAGS = -std=c++11 -g -I. -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -DVERBOSE
  LINKER = g++ -std=c++11 -o
else
  CFLAGS = -std=c++11 -I. -O3
  LINKER = g++ -std=c++11 -O3 -o
endif

LFLAGS = -I.

# directories
SRCDIR = src
OBJDIR = obj
BINDIR = bin

# our program
TARGET = cafexp

SOURCES := $(wildcard $(SRCDIR)/*.cpp)
#INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
rm = rm -f

$(BINDIR)/$(TARGET): $(OBJECTS) $(OBJDIR)/main.o
	$(LINKER) $@ $(LFLAGS) $(OBJECTS) $(OBJDIR)/main.o
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

$(OBJDIR)/main.o: main.cpp
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

$(OBJDIR)/test.o: test.cpp
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

.PHONY: clean
clean:
	$(rm) $(OBJECTS) $(OBJDIR)/main.o
	@echo "Cleanup complete!"

.PHONY: remove
remove: clean
	@$(rm) $(BINDIR)/$(TARGET)
	@echo "Executable removed!"

TESTDIR = $(BINDIR)
TEST_TARGET = runtests
TESTLIBS=-lCppUTest -lCppUTestExt
$(TESTDIR)/$(TEST_TARGET): $(OBJECTS) $(OBJDIR)/test.o
	$(LINKER) $@ $(LFLAGS) $(OBJECTS) $(OBJDIR)/test.o $(TESTLIBS)
	@echo "Test linking complete!"

.PHONY: test
test: $(TESTDIR)/$(TEST_TARGET)


# -- Notes for later, when I forget GNU make syntax: -- #

# := means a simply expanded variable, which behave more or
# less like a normal programming language variable.
# References to other variables in a simply expanded variable
# are expanded just once at the moment of variable creation,
# and that will be its final value.

# Wildcard expansion happens automatically (w/o the need to
# invoke a function for it) in rules, but not when variables
# are set. In those cases, one needs to call the 'wildcard' function.

# The symbol '%' in wildcard expansions means any character any
# number of times. And if both the searched and replacement patterns
# have '%', the '%' in replacement will match the '%' in the search
# pattern.

# The form '$(var:a=b)' or '${var:a=b}' is called a substitution
# reference. It is an abbreviation of the 'patsubst' function, and
# replaces every 'a' at the end of the values of 'var' (has to be at
# the end) with 'b'. Another way to use substitution references is
# with '$(var:%.o=%.cpp)', for instance; here, this expression will
# evaluate to the names of all .o files, but replacing .o with .cpp.

# The form 'targets: target-pattern: prereq-pattern' is called a
# static pattern rule. This is more general than a rule with multiple
# targets that must share the same prerequisites. Here, the prereqs
# are named according to the target names. If we do not use this,
# and use something like $(OBJECTS): src/*.cpp, then we are saying
# each and all .o files will require all each and all .cpp files,
# which will lead to errors.

# -I. means that g++ should look for .h files in the same folder where
# the .cpp files are.

# The '@' symbol at the start of a line suppresses the echoing of
# the line.

# The -c flag tells g++ to compile object files.

# The -o flag tells g++ what the object files should be named.

# $< is replaced by the first item in the list of prerequisites (which
# starts after the ':' in a rule).

# $@ is replaced by the rule's target name.

# .PHONY specifies a phony target. This is good practice to avoid
# conflict with a file with the same name as the target, say, a file
# called 'clean'. In such case, the 'clean' file would always be
# up-to-date and the recipe of 'clean' would never be executed.
# And if a file named 'clean' is never created, then a rule such as
# clean: blabla (with some recipe) will ALWAYS be executed.

# ---------------------------------------------------- #
