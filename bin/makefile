# Define source code directories
SOURCE_DIR = ../src
WTP_DIR = $(SOURCE_DIR)/wtp
WTP_OPTIMIZE_DIR = $(SOURCE_DIR)/wtp_optimize
AUTODOSE_DIR = $(SOURCE_DIR)/autodose
BORG_DIR = $(SOURCE_DIR)/borg

SOURCES =                              \
	$(SOURCE_DIR)/main.cpp               \
	$(WTP_DIR)/adj_fact.cpp              \
	$(WTP_DIR)/alm_rmv2.cpp              \
	$(WTP_DIR)/basin3.cpp                \
	$(WTP_DIR)/breakpt.cpp               \
	$(WTP_DIR)/chemical.cpp              \
	$(WTP_DIR)/cl2decy3.cpp              \
	$(WTP_DIR)/clo2.cpp                  \
	$(WTP_DIR)/coagdbp.cpp               \
	$(WTP_DIR)/coef.cpp                  \
	$(WTP_DIR)/ct.cpp                    \
	$(WTP_DIR)/data_wtp.cpp              \
	$(WTP_DIR)/dbpslect.cpp              \
	$(WTP_DIR)/distdbp3.cpp              \
	$(WTP_DIR)/ec.cpp                    \
	$(WTP_DIR)/fcl_rmv2.cpp              \
	$(WTP_DIR)/fdoublef.cpp              \
	$(WTP_DIR)/filter3.cpp               \
	$(WTP_DIR)/ftab.cpp                  \
	$(WTP_DIR)/gacmbdbp.cpp              \
	$(WTP_DIR)/gac_rmv.cpp               \
	$(WTP_DIR)/globals.cpp               \
	$(WTP_DIR)/influent.cpp              \
	$(WTP_DIR)/list_wtp.cpp              \
	$(WTP_DIR)/mdrw1dbp.cpp              \
	$(WTP_DIR)/mdrw2dbp.cpp              \
	$(WTP_DIR)/mfuf_rmv.cpp              \
	$(WTP_DIR)/ncstr.cpp                 \
	$(WTP_DIR)/nf_rmv.cpp                \
	$(WTP_DIR)/open_wtp.cpp              \
	$(WTP_DIR)/owdbp.cpp                 \
	$(WTP_DIR)/ozone.cpp                 \
	$(WTP_DIR)/phchange.cpp              \
	$(WTP_DIR)/readline.cpp              \
	$(WTP_DIR)/read_wtp.cpp              \
	$(WTP_DIR)/res_time.cpp              \
	$(WTP_DIR)/runmodel.cpp              \
	$(WTP_DIR)/run_wtp.cpp               \
	$(WTP_DIR)/rwdbp.cpp                 \
	$(WTP_DIR)/save_wtp.cpp              \
	$(WTP_DIR)/soft_rmv.cpp              \
	$(WTP_DIR)/solids.cpp                \
	$(WTP_DIR)/strings.cpp               \
	$(WTP_DIR)/struct.cpp                \
	$(WTP_DIR)/thm_wtp.cpp               \
	$(WTP_DIR)/uptable.cpp               \
	$(WTP_DIR)/writewtp.cpp              \
	$(WTP_OPTIMIZE_DIR)/read_csv.cpp \
	$(WTP_OPTIMIZE_DIR)/wtp_optimize.cpp    \
	$(WTP_OPTIMIZE_DIR)/wtp_problem.cpp     \
	$(BORG_DIR)/borg.cpp                 \
	$(BORG_DIR)/mt19937ar.cpp            \
	$(AUTODOSE_DIR)/auto_dose.cpp          \
	$(AUTODOSE_DIR)/extrema.cpp               \
	$(AUTODOSE_DIR)/rootfind_and_mod_dose.cpp \


OBJECTS=$(SOURCES:.cpp=.o)  
# CPP=g++
CPP=g++ -std=c++11  # use C++11
# CPPFLAGS=-g -c -O3 -Wall -I. -I$(SOURCE_DIR) -I$(AUTODOSE_DIR) -I$(BORG_DIR) -I$(WTP_DIR) -I$(WTP_OPTIMIZE_DIR)  # for best performance
CPPFLAGS=-g -c -O0 -Wall -I. -I$(SOURCE_DIR) -I$(AUTODOSE_DIR) -I$(BORG_DIR) -I$(WTP_DIR) -I$(WTP_OPTIMIZE_DIR)  # for debugging

LIBS=-lm
EXECUTABLE=wtp-optimize.exe

all: $(SOURCES) $(EXECUTABLE)
	rm -f $(SOURCE_DIR)/*.o
	rm -f $(WTP_DIR)/*.o 
	rm -f $(WTP_OPTIMIZE_DIR)/*.o 
	rm -f $(AUTODOSE_DIR)/*.o 
	rm -f $(BORG_DIR)/*.o 

$(EXECUTABLE): $(OBJECTS)
	$(CPP) $(OBJECTS) -o $@ $(LIBS)  

.cpp.o: 
	$(CPP) $(CPPFLAGS) $^ -o $@  

clean:
	rm -f $(SOURCE_DIR)/*.o $(EXECUTABLE)
	rm -f $(WTP_DIR)/*.o 
	rm -f $(WTP_OPTIMIZE_DIR)/*.o 
	rm -f $(AUTODOSE_DIR)/*.o 
	rm -f $(BORG_DIR)/*.o 

