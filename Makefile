PROJ_NAME = fem
CC = gcc

OBJ_DIR = obj
OBJ = $(SRC:src/%.c=obj/%.o)
SRC = $(wildcard src/*.c)
LINK = -lgsl -lgslcblas -lm
CFLAGS = -Wall
INCLUDE = include


$(PROJ_NAME): $(OBJ_DIR) $(OBJ)
	$(CC) -o $@ $(OBJ)  $(LINK)

$(OBJ_DIR)/main.o: ./src/main.c
	$(CC) -o $@ -c $< $(CFLAGS) -I$(INCLUDE)

$(OBJ_DIR)/%.o : ./src/%.c ./include/%.h
	$(CC) -o $@ -c $< $(CFLAGS) -I$(INCLUDE)  $(LINK)

.PHONY: clean
clean:
	rm -f ./$(OBJ_DIR)/*.o $(PROJ_NAME)

$(OBJ_DIR): 
	mkdir $(OBJ_DIR)

nosx = 51
nosy = 51
limInfx = 0.0
limSupx = 1.0
limInfy = 0.0
limSupy = 1.0
erro = 1e-6
iteracoes = 10000

parametros = $(nosx) $(nosy) $(limInfx) $(limSupx) $(limInfy) $(limSupy) $(erro) $(iteracoes)
run: $(PROJ_NAME)
	rm -rf resultados/*
	./$(PROJ_NAME) $(parametros)
#	make imagem > /dev/null
#	make plot > /dev/null

imagem:
	gnuplot grafico.plt

plot:
	gnuplot plot.plt 2> /dev/null

