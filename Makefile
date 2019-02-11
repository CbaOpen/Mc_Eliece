PROJECT = MC_ELIECE_Project
EXE_NAME = mceliece

SRC_PATH = src/
INC_PATH = inc/
OBJ_PATH = obj/
EXE_PATH = bin/

EXEC = $(EXE_PATH)$(EXE_NAME)
SRC = $(shell find $(SRC_PATH)*.c)
INC = $(shell find $(INC_PATH)*.h)
OBJ = $(SRC:$(SRC_PATH)%.c=$(OBJ_PATH)%.o)

## Compilation ................................................................:

CC = gcc
INC_FLAGS    = -I$(INC_PATH)
PROJ_CFLAG   = -c -g3 -Wall
PROJ_LDFLAGS = -Wall #-lgmp -lm

CFLAGS  = $(PROJ_CFLAG) $(INC_FLAGS)
LDFLAGS = $(PROJ_LDFLAG)

## Lancement ..................................................................:

ARGS = -h

# Cibles =======================================================================

.PHONY : clean mrproper

## Lancement ..................................................................:

run : compil 
	@echo "--> Lancement de '$(EXEC)' :"
	$(EXEC) $(ARGS)

## Compilation ................................................................:

compil : $(EXEC)

$(EXEC) : $(OBJ) 
	@echo "--> Édition des liens dans '$(EXEC)' :"
	$(CC) $^ -o $(EXEC) $(LDFLAGS)

$(OBJ_PATH)%.o : $(SRC_PATH)%.c Makefile
	@echo "--> Compilation de '$<' :"
	$(CC) -c $< -o $@ $(CFLAGS)

## Dépendances ................................................................:

-include $(OBJ:%.o=%.d)

## Nettoyage ..................................................................:

clean :
	@echo "--> Suppression des fichier temporaires de $(PROJECT) :"
	rm -f $(OBJ_PATH)*.o $(OBJ_PATH)*.d $(SRC_PATH)*~ $(INC_PATH)*~

mrproper: clean
	@echo "--> Suppression de l'executable de $(PROJECT) :"
	rm -f $(EXEC)

## Debugger & Profiler ........................................................:

gdb : compil
	@echo "--> Debbugage avec $@ :"
	$@ -q --args $(EXEC) $(ARGS)

valgrind-p1 : compil
	@echo "--> Debbugage avec $@ (profile 1) :"
	valgrind --tool=memcheck --leak-resolution=high \
	    --show-possibly-lost=yes --show-reachable=yes $(EXEC) $(ARGS)

valgrind-p2 : compil
	@echo "--> Debbugage avec $@ (profile 2) :"
	valgrind --tool=memcheck --leak-resolution=high --leak-check=full \
	    --show-possibly-lost=yes --show-reachable=yes --track-origins=yes \
	    $(EXEC) $(ARGS)
