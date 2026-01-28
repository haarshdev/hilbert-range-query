# Copyright (C) Jonathan Lawder 2001-2011

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		EXECUTABLES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARGET1      =  a.exe
TARGET2      =  b.exe
DEMO	     =  demo.exe
TARGETj      =  j.exe

#..............................................................................
#		IF FDL NOT ENABLED			FDLFDLFDLFDLFDL!!!!!!!!
#..............................................................................
#All: $(TARGET1) $(TARGET2)
All: $(DEMO)# $(TARGET2)

FDL          =  0
LTC          =  0
NON_KEY_INFO =	0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		DIRECTORIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ROOT_DIR	=	..\\
B_DIR	=	$(ROOT_DIR)btree\\
D_DIR	=	$(ROOT_DIR)db\\
H_DIR	=	$(ROOT_DIR)hilbert\\
L_DIR	=	$(ROOT_DIR)ltc\\
U_DIR	=	$(ROOT_DIR)utils\\
T_DIR	=	$(ROOT_DIR)tests\\

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		DEFAULT #defines: may be over-ridden by environment variables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# these are spare macros
EXTRA1       =	DEV
EXTRA2       =	NONE2
EXTRA3       =	NONE3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		DEFINES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DEFINES      =	-D$(EXTRA1) -D$(EXTRA2) -D$(EXTRA3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		COMPILATION FLAGS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GPROF        =	#-xpg
DBX          =  -Wall #-g
OPTIMIZE     =  #-fast -native

FLAGS        =  $(GPROF) $(DBX) $(OPTIMIZE) # -xarch=v9a (64 bit)

C_FLAGS      =	-v -I$(I_DIR) $(FLAGS)
T_FLAGS      =	$(C_FLAGS) $(DEFINES) -o	# target flags
O_FLAGS      =	$(C_FLAGS) $(DEFINES) -c	# object flags

COMPILER     =  gxx	# or gcc?
COMPILER2    =  gxx

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		OBJECT DEFINITIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OBJECTS1	=	btree.o utils.o test5.o
OBJECTS2	=	btree.o db.o buffer.o page.o hilbert.o utils.o test2.o
DEMO_OBJ	=	btree.o db.o buffer.o page.o query.o hilbert.o utils.o demo.o
OBJECTSj	=	db.o buffer.o page.o utils.o testj.o

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		TARGET DEFINITIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

$(TARGET1):	$(OBJECTS1)
		$(COMPILER2) $(T_FLAGS) $(TARGET1) $(OBJECTS1)

$(TARGET2):	$(OBJECTS2)
		$(COMPILER2) $(T_FLAGS) $(TARGET2) $(OBJECTS2)

$(DEMO):	$(DEMO_OBJ)
		$(COMPILER2) $(T_FLAGS) $(DEMO) $(DEMO_OBJ)

$(TARGETj):	$(OBJECTSj)
		$(COMPILER2) $(T_FLAGS) $(TARGETj) $(OBJECTSj)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		DEPENDENCIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test5.o:	$(ROOT_DIR)gendefs.h $(B_DIR)btree.h $(U_DIR)utils.h \
		$(T_DIR)test5.cc
		$(COMPILER) $(O_FLAGS) $(T_DIR)test5.cc

test2.o:	$(ROOT_DIR)gendefs.h $(B_DIR)btree.h $(U_DIR)utils.h \
		$(T_DIR)test2.cc
		$(COMPILER) $(O_FLAGS) $(T_DIR)test2.cc

demo.o:		$(ROOT_DIR)gendefs.h $(B_DIR)btree.h $(U_DIR)utils.h \
		$(D_DIR)db.h $(D_DIR)buffer.h $(D_DIR)page.h\
		$(T_DIR)demo.cc
		$(COMPILER) $(O_FLAGS) $(T_DIR)demo.cc

testj.o:	$(ROOT_DIR)gendefs.h $(U_DIR)utils.h \
		$(D_DIR)db.h buffer.h page.h \
		$(T_DIR)testj.cc
		$(COMPILER) $(O_FLAGS) $(T_DIR)testj.cc

btree.o:	$(ROOT_DIR)gendefs.h $(B_DIR)btree.h $(U_DIR)utils.h \
		$(B_DIR)btree.cc
		$(COMPILER) $(O_FLAGS) $(B_DIR)btree.cc

db.o:		$(ROOT_DIR)gendefs.h $(B_DIR)btree.h $(D_DIR)db.h $(D_DIR)db.cc
		$(COMPILER) $(O_FLAGS) $(D_DIR)db.cc

buffer.o:	$(ROOT_DIR)gendefs.h $(D_DIR)db.h $(D_DIR)buffer.h \
		$(D_DIR)page.h $(D_DIR)buffer.cc
		$(COMPILER) $(O_FLAGS) $(D_DIR)buffer.cc

page.o:		$(ROOT_DIR)gendefs.h $(D_DIR)db.h $(D_DIR)page.h \
		$(D_DIR)page.cc
		$(COMPILER) $(O_FLAGS) $(D_DIR)page.cc

query.o:	$(ROOT_DIR)gendefs.h $(D_DIR)db.h $(D_DIR)query.cc
		$(COMPILER) $(O_FLAGS) $(D_DIR)query.cc

hilbert.o:	$(ROOT_DIR)gendefs.h $(U_DIR)utils.h $(H_DIR)hilbert.cc
		$(COMPILER) $(O_FLAGS) $(H_DIR)hilbert.cc

utils.o:	$(ROOT_DIR)gendefs.h $(U_DIR)utils.cc
		$(COMPILER) $(O_FLAGS) $(U_DIR)utils.cc