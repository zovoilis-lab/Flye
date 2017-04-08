PREFIX		=	.
CXXFLAGS	= 	-Wall -W -pedantic-errors -Wno-unused-parameter # -Werror
CXXFLAGS	+=  -Wmissing-braces -Wparentheses -ansi

# Directories

PILDIR		=	$(PREFIX)/pilal
SIMDIR		=	$(PREFIX)/simplex

CPPFLAGS    =   -I$(PILDIR)/include -I$(SIMDIR)/include

# Cleaning options
RM			= 	rm -rf

# Files
SRCF		=	$(wildcard $(SRCDIR)/*.cc) $(wildcard $(PILDIR)/src/*.cc) $(wildcard $(SIMDIR)/src/*.cc)
DEPS		=	$(SRCF:.cc=.d*)
OBJS		=	$(SRCF:.cc=.o)
LIB			=   libcpplex.so
#BINS		=	$(BINDIR)/solver

# Targets
# all library clean wipe
.PHONY:	all clean wipe

all: CXXFLAGS += -O3 -DNDEBUG
all: libcpplex

debug: CXXFLAGS += -Og -g -D_GLIBCXX_DEBUG
debug: libcpplex

#$(CXX) $^ -o $@
libcpplex: $(OBJS)
	ar rs $(LIB) $(OBJS)

clean:
	$(RM) $(OBJS) $(DEPS)

wipe:
	$(RM) $(OBJS) $(DEPS) $(BINS)

%.d: %.cc
	@set -e; rm -f $@; \
	$(CXX) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$ \

include $(SRCF:.cc=.d)
