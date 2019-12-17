
#include<stdio.h>
#include <string.h>
//===================================================================

void scriptFigGnuplot(int nosx, int nosy, double limInfx, double limSupx, double limInfy, double limSupy, int n){
	FILE *h;
	char str[18]="figSolucao";
	sprintf(str, "figSolucao%03d.plt",n);
	h = fopen(str,"w");

	fprintf(h,"%s\n","reset");
	fprintf(h,"set xrange [%.2lf:  %.2lf]\n",limInfx,limSupx);
	fprintf(h,"set yrange [%.2lf:  %.2lf]\n",limInfy,limSupy);
//	fprintf(h,"%s\n","set zrange [0.0:  3.5]");// f variavel b variavel
//	fprintf(h,"%s\n","set xtics 0.2");
//	fprintf(h,"%s\n","set ytics 0.2");
//	fprintf(h,"%s\n","set style data lines");
	//fprintf(h,"%s\n","set terminal postscript");
	fprintf(h,"%s\n","set key left top");
	fprintf(h,"%s\n","set xlabel 'x'");
	fprintf(h,"%s\n","set ylabel 'y'");
	fprintf(h,"%s\n","set pm3d map");

//	fprintf(h,"%s\n","set cbrange [0:50]"); //fonte constante
//	fprintf(h,"%s\n","set cbrange [0:3.5]"); //f variavel b variavel

//	fprintf(h,"%s %d, %d, %d\n","set dgrid3d",nosy,nosx,2);
	fprintf(h,"%s\n","set term png");
	fprintf(h,"set output 'figura%03d.png'\n",n);
	fprintf(h,"splot \"galerkin%d.txt\" notitle\n",n);
	fclose(h);
}
//===================================================================
void scriptLatex(int n)
{
	FILE *f;
	int i;
	f = fopen("join.tex","w");
	fprintf(f, "\\documentclass{article}\n\n");
	fprintf(f, "\\usepackage{graphicx}\n");
	fprintf(f, "\\begin{document}\n\n");
	for (i=0; i<n; i++)
	{
		fprintf(f, "\\centering\n");
		fprintf(f, "\\includegraphics[width=0.9\\textwidth]{figura%d.eps}\n", i);
		fprintf(f, "\\newpage\n\n");
	}
	fprintf(f, "\\centering\n");
	fprintf(f, "\\includegraphics[width=0.9\\textwidth]{figura%d.eps}\n", n);
	fprintf(f, "\\end{document}\n");
	fclose(f);
}
