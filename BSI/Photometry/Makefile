
#@author: Conor O'Toole
#		 University College Dublin,
#		 Ireland
#		 
#@supervisor: Dr. Morgan Fraser
#			 University College Dublin,
#			 Ireland
#			 
#			 
#		 Supported by University College Dublin Seed Funding Scheme 2017


# Put the directory above 'DNest4' into an environment variable
# called DNEST4_PATH for this to work.

CC = g++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -pedantic -DNDEBUG
LIBS = -ldnest4 -lpthread
DNEST4_PATH = /c/Program Files/

default:
	$(CC) -I$(DNEST4_PATH) $(CXXFLAGS) -c *.cpp
	$(CC) -pthread -L$(DNEST4_PATH)/DNest4/code -o main *.o $(LIBS)
	rm *.o

