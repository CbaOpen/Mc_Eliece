PROJECT = MC_ELIECE_Cryptosystem
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
PROJ_CFLAG   = -c -g3 -Wall -lflint -lmpfr -lgmp -lpthread -lm
PROJ_LDFLAGS = -Wall -lflint -lmpfr -lgmp -lpthread -lm

CFLAGS  = $(INC_FLAGS) $(PROJ_CFLAG)
LDFLAGS = $(PROJ_LDFLAGS)

## Lancement ..................................................................:

ARGS =

# Cibles =======================================================================

.PHONY : clean mrproper

## Compilation ................................................................:

compil : $(EXEC)

$(EXEC) : $(OBJ) 
	@echo "--> Édition des liens dans '$(EXEC)' :"
	$(CC) $^ -o $(EXEC) $(LDFLAGS)

$(OBJ_PATH)%.o : $(SRC_PATH)%.c Makefile
	@echo "--> Compilation de '$<' :"
	$(CC) -c $< -o $@ $(CFLAGS)

## Lancement ..................................................................:

run : compil 
	@echo "--> Lancement de '$(EXEC)' :"
	$(EXEC) $(ARGS)

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
