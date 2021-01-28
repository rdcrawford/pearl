
# Get the files to compile into pearl
target = pearl
src :=  BlastAlignment.cpp  BlastData.cpp  BlastResults.cpp  InputParser.cpp \
  Subject.cpp  GenomeData.cpp Genome.cpp BioSeq.cpp pearl.cpp
objects:= $(addsuffix .o, $(basename $(notdir $(src))))

libs = -lstdc++fs
cxx := g++
flags = -Wall -std=c++17

# compile object files:
%.o: %.cpp
	$(cxx) $(includes) $(flags) -c $< -o $@ $(libs)
	@echo $@

all: $(target)

# build exe
$(target): $(objects)
	$(cxx) $(flags) $(objects) -o $(target) $(libs)

clean:
	rm -f $(objects) $(target)

cleanexe:
	rm -f $(target)

print:
	@echo $(objects)
