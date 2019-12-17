reset
set pm3d map at s
set xlabel 'x'
set ylabel 'y'
set view 66,342,1,1
unset key
#set view 90, 270,1,1
set term png

t = system("ls ./resultados/dados | grep -c resultado")
system "mkdir -p ./resultados/imagens"
do for [i=0:t-1] {
	name = sprintf('iteração %03d', i)
	set title name
  outfile = sprintf('./resultados/imagens/figura%.0f.png',i)
  set output outfile
  inputfile = sprintf('./resultados/dados/resultado%.0f.txt',i)
  splot inputfile
}

