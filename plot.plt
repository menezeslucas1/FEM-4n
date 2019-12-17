reset

#set terminal qt 0 position 100, 50
#set style data lines
#set xlabel 'iteração'
#set ylabel 'log Ajuste'
#unset key
#inputfile = sprintf('./resultados/ajuste.txt')
#plot inputfile #title "data"

#set terminal qt 1 position 742, 50
set pm3d map at s
set xlabel 'x'
set ylabel 'y'
unset key
set view 66,342,1,1
t = system("ls ./resultados/dados | grep -c resultado")
inputfile = sprintf('./resultados/dados/resultado%.0f.txt',t-1)
#t = system("ls ./resultados/dados")
#inputfile = sprintf('./resultados/dados/%s',t)
splot inputfile

pause -1

