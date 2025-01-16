target = Gt_Gw

SOURCES = correlator_mode0.cc correlator_mode1.cc correlator_Gw.cc main.cc	
OBJECTS = $(SOURCES:.cc=.o) 

CXXFLAGS=-O3 

default: $(target) 

$(target): $(OBJECTS) 
	$(CXX) -o $@ $^

clean:
	-rm -f *.o $(target) 

%.o: %.cc 
	$(CXX) $(CXXFLAGS.o)  -c $< -o $@

depend:
	makedepend -Y. -o.o -- $(SOURCES) 
# DO NOT DELETE

