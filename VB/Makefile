TARGET	= VB

CC		= g++
CFLAGS	= -W -Wall
LINKER	= g++
LFLAGS	= -lm -W -Wall -o 
OBJDIR	= obj
SRCDIR	= src
BINDIR	= bin

SOURCES := $(wildcard $(SRCDIR)/*.cpp)
INCLOUDS := $(wildcard $(SRCDIR)/*.h)
OBJECTS	:= $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

all: $(BINDIR)/$(TARGET)

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

$(BINDIR)/$(TARGET): $(BINDIR) $(OBJECTS)
	$(LINKER) $(LFLAGS) $@ $(OBJECTS)
	@echo "Linking complete!"

$(OBJDIR):
	mkdir $(OBJDIR)

$(BINDIR):
	mkdir $(BINDIR)

clean:
	rm -rf $(OBJDIR) $(BINDIR)/$(TARGET)
