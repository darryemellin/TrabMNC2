#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
Autores:
    Darrye Roberto da Silva Mellin
    Willyan Pizarro Pretes
    Bernardo Betete da Silva Fonseca
    */

void iniciaTabela(int numeroPontos,double tabela[numeroPontos][numeroPontos]){
    int coluna=0,linha=0;
    for(coluna=0;coluna<numeroPontos;coluna++){
        for(linha=0;linha<numeroPontos;linha++){
            tabela[coluna][linha]=0;
        }
    }
}

void Newton(int numeroPontos,double tabelaPontos[2][numeroPontos],double pontoDesejado){
    double tabelaAuxiliar[numeroPontos][numeroPontos];
    iniciaTabela(numeroPontos,tabelaAuxiliar);
    int coluna=0,linha=0;
    int contadorLinha=numeroPontos;
    int pulox=0;
    double resultado=0,auxiliar=1;
    for(coluna=0;coluna<numeroPontos;coluna++){
        for(linha=0;linha<contadorLinha;linha++){
            if(coluna==0){
                tabelaAuxiliar[0][linha]=tabelaPontos[1][linha];
            }
            else{
                tabelaAuxiliar[coluna][linha]= (tabelaAuxiliar[coluna-1][linha+1]-tabelaAuxiliar[coluna-1][linha])/(tabelaPontos[0][linha+pulox]-tabelaPontos[0][linha]);
            }
        }
        contadorLinha--;
        pulox++;
    }
    //y0 + (x-x0)*y1 + (x-x0)(x-x1)*y2 + (x-x0)(x-x1)(x-x3)*y3 ...
    contadorLinha = numeroPontos-1;
    linha=0;
    coluna=0;
    for(coluna=0;coluna<numeroPontos-1;coluna++){//repeticoes
        for(linha=contadorLinha-1;linha>=0;linha--){// de tras pra frente
            auxiliar *=(pontoDesejado-tabelaPontos[0][linha]);
        }
        resultado += auxiliar*tabelaAuxiliar[contadorLinha][0];//y respectivo
        auxiliar=1;
        contadorLinha--;
    }
    resultado += tabelaAuxiliar[0][0];

    printf("Resultado no ponto %.4lf: %.4lf",pontoDesejado,resultado);
}

int fatorial(int valor){
    int resposta=1;
    int i;
    for(i=valor;i>0;i--){
        resposta *= i;
    }
    return resposta;
}

void NewtonGregory(int numeroPontos,double tabelaPontos[2][numeroPontos],double pontoDesejado){
    double tabelaAuxiliar[numeroPontos][numeroPontos];
    double distancia=0;
    distancia = tabelaPontos[0][1]-tabelaPontos[0][0];
    iniciaTabela(numeroPontos,tabelaAuxiliar);
    int coluna=0,linha=0;
    int contadorLinha=numeroPontos;
    int pulox=0;
    double resultado=0,auxiliar=1;
    for(coluna=0;coluna<numeroPontos;coluna++){
        for(linha=0;linha<contadorLinha;linha++){
            if(coluna==0){
                tabelaAuxiliar[0][linha]=tabelaPontos[1][linha];
            }
            else{
                tabelaAuxiliar[coluna][linha]= (tabelaAuxiliar[coluna-1][linha+1]-tabelaAuxiliar[coluna-1][linha]);
            }
        }
        contadorLinha--;
        pulox++;
    }
    //y0 + (x-x0)*y1/h + (x-x0)(x-x1)*y2/2!h^2 + (x-x0)(x-x1)(x-x3)*y3/3!^3 ...
    contadorLinha = numeroPontos-1;
    linha=0;
    coluna=0;
    for(coluna=0;coluna<numeroPontos-1;coluna++){//repeticoes
        for(linha=contadorLinha-1;linha>=0;linha--){// de tras pra frente
            auxiliar *=(pontoDesejado-tabelaPontos[0][linha]);
        }
        resultado += auxiliar*(tabelaAuxiliar[contadorLinha][0]/(fatorial(contadorLinha)*pow(distancia,contadorLinha)));//y respectivo com a conta
        auxiliar=1;
        contadorLinha--;
    }
    resultado += tabelaAuxiliar[0][0];

    printf("Resultado no ponto %.4lf: %.4lf",pontoDesejado,resultado);
}

double CoefDeterminacao(int numeroPontos,double tabelaPontos[2][numeroPontos],double yajustado[numeroPontos]){
    double erroquad=0,erro=0;
    double somay=0, somayquad=0,resposta=0;
    int i=0;
    for(i=0;i<numeroPontos;i++){
        erro = tabelaPontos[1][i]-yajustado[i];
        erroquad += erro*erro;
        somay += tabelaPontos[1][i];
        somayquad+= tabelaPontos[1][i]*tabelaPontos[1][i];
    }
    
    resposta = ((numeroPontos*erroquad)/(numeroPontos*somayquad - somay*somay));
    return 1-resposta;
    
}

void AjusteReta(int numeroPontos,double tabelaPontos[2][numeroPontos],double a0,double a1,double yajustado[numeroPontos],double coefDeterminacao){
    double somax=0,somay=0,somaxy=0,somaxquad=0;
    int i=0;
    for(i=0;i<numeroPontos;i++){
        somax += tabelaPontos[0][i];
        somay += tabelaPontos[1][i];
        somaxy += tabelaPontos[0][i]*tabelaPontos[1][i];
        somaxquad += tabelaPontos[0][i]*tabelaPontos[0][i];
    }
    a1 = (numeroPontos*somaxy-somax*somay)/(numeroPontos*somaxquad-(somax*somax));
    a0 = (somay-a1*somax)/numeroPontos;

    for(i=0;i<numeroPontos;i++){
        yajustado[i]= a0 + a1*tabelaPontos[0][i];
    }

    coefDeterminacao = CoefDeterminacao(numeroPontos,tabelaPontos,yajustado);

    printf("equacao\ny = %.4lf + (%.4lf)*x\n",a0,a1);
    printf("coeficiente de determinacao: %.4lf\n",coefDeterminacao);


}

void RotinaDecomposicaoLU (int ordemMatriz, double matrizA[ordemMatriz][ordemMatriz], double termosIndependentes[ordemMatriz], double vetorSolucao[ordemMatriz]){
	double det;
    int a, b, i, j, k,aux;
	double soma = 0, lsoma = 0, lsub = 0, usoma = 0, usub = 0;
	double matrizu[ordemMatriz][ordemMatriz], matrizl[ordemMatriz][ordemMatriz], vetaux[ordemMatriz];

	for (i = 0; i < ordemMatriz; i++)
    {
		for (j = 0; j < ordemMatriz; j++)
        {
            matrizu[i][j] = 0;
            matrizl[i][j] = 0;
            vetaux[i] = 0;
            if (i == j){
                matrizl[i][j] = 1;
            }
	   }
	}

	for (a = 0; a < ordemMatriz; a++)
    {
		for (b = 0; b < ordemMatriz; b++)
        {
            if(a <= b)
            {
                usoma = 0;
                for (k = 0; k < a; k++)
                    usoma = usoma + (matrizl[a][k] * matrizu[k][b]);

                matrizu[a][b] = matrizA[a][b] - usoma;

            }
            if(a > b)
            {
                usub = 0;
                for (k = 0; k < b; k++)
                    usub = usub + (matrizl[a][k] * matrizu[k][b]);

                matrizl[a][b] = (matrizA[a][b] - usub)/matrizu[b][b];

            }
        }
	}
	
    for (a = 0; a < ordemMatriz; a++)
    {
        soma = 0;	
        for (b = 0; b < a; b++){
            soma = soma - matrizl[a][b] * vetaux[b];
        }
        if (soma == 0) vetaux[a] = termosIndependentes[a]/matrizl[a][a];
        else vetaux[a] = (termosIndependentes[a] + soma)/matrizl[a][a];
    }  
	  
    for (a = ordemMatriz-1; a >= 0; a--)
    {
        soma = 0;	
        for (b = ordemMatriz-1; b > a; b--){
            soma = soma - matrizu[a][b] * vetorSolucao[b];
        }
        if (soma == 0) vetorSolucao[a] = vetaux[a]/matrizu[a][a];
        else vetorSolucao[a] = (vetaux[a] + soma)/matrizu[a][a];
	}  
  
}

void AjustePolinomial(int numeroPontos,int grau,double tabelaPontos[2][numeroPontos],double coeficientes[],double yajustado[numeroPontos],double coefDeterminacao){
    double tabelaAuxiliar[grau+1][grau+1];
    double termosIndependentes[grau+1];
    int linha,coluna=1,i,elevacao=0;
    iniciaTabela(grau+1,tabelaAuxiliar);
    tabelaAuxiliar[0][0]=numeroPontos;

    for(linha=0;linha<=grau;linha++){
        for(coluna=0;coluna<=grau;coluna++){
            for(i=0;i<numeroPontos;i++){
                if(linha==0 && coluna==0){
                    tabelaAuxiliar[0][0]=numeroPontos;
                }
                else{
                    tabelaAuxiliar[linha][coluna] += pow(tabelaPontos[0][i],elevacao+coluna);
                }
            }
        }
        elevacao++;
    }

    for(i=0;i<grau+1;i++){
        termosIndependentes[i]=0;
        for(linha=0;linha<numeroPontos;linha++){
            termosIndependentes[i] += tabelaPontos[1][linha]*pow(tabelaPontos[0][linha],i); 
        }
    }

    RotinaDecomposicaoLU(grau+1,tabelaAuxiliar,termosIndependentes,coeficientes);

    double auxsoma=0,potencia=0;
    for(i=0;i<numeroPontos;i++){
        for(linha=0;linha<grau+1;linha++){
            potencia =pow(tabelaPontos[0][i],linha);
            auxsoma += coeficientes[linha]*potencia;
        }
        yajustado[i] = auxsoma;
        auxsoma=0;
    }

    coefDeterminacao = CoefDeterminacao(numeroPontos,tabelaPontos,yajustado);
    printf("coeficientes\n");
    for(i=0;i<grau+1;i++){
        printf("a%d = %.4lf\n",i,coeficientes[i]);
    }
    printf("coeficiente de determinacao: %.4lf\n",coefDeterminacao);

}

void AjusteExponencial(int numeroPontos,double tabelaPontos[2][numeroPontos],double a,double b,double yajustado[numeroPontos],double coefDeterminacao){
    double ymodificado[numeroPontos];
    double a1=0,a0=0;
    double somax=0,somay=0,somaxy=0,somaxquad=0;
    int i=0;
    for(i=0;i<numeroPontos;i++){
        ymodificado[i] = log(tabelaPontos[1][i]);
    }

    for(i=0;i<numeroPontos;i++){
        somax += tabelaPontos[0][i];
        somay += ymodificado[i];
        somaxy += tabelaPontos[0][i]*ymodificado[i];
        somaxquad += tabelaPontos[0][i]*tabelaPontos[0][i];
    }
    a1 = (numeroPontos*somaxy-somax*somay)/(numeroPontos*somaxquad-(somax*somax));
    a0 = (somay-a1*somax)/numeroPontos;

    for(i=0;i<numeroPontos;i++){
        yajustado[i]= a0 + a1*tabelaPontos[0][i];
    }

    for(i=0;i<numeroPontos;i++){
        tabelaPontos[1][i]=ymodificado[i];
    }

    coefDeterminacao = CoefDeterminacao(numeroPontos,tabelaPontos,yajustado);

    a= exp(a0);
    b= exp(a1);

    printf("equacao reta\ny = %.4lf + (%.4lf)*x\n",a0,a1);
    printf("equacao exponencial\ny = %.4lf*(%.4lf^x)\n",a,b);
    printf("coeficiente de determinacao: %.4lf\n",coefDeterminacao);
}


int main(int argc, char** argv) {
    int numeroPontosTabelados;
    printf("Informe a quantidade de pontos tabelados\n");
    scanf("%d",&numeroPontosTabelados);

    double tabelaPontos[2][numeroPontosTabelados];
    double a1,a0,a,b,yajustado[numeroPontosTabelados],coefDeterminacao;
    int i;
    printf("Popule a parte X\n");
    for(i=0;i<numeroPontosTabelados;i++){
        printf("x%d: ",i);
        scanf("%lf",&tabelaPontos[0][i]);
    }

    printf("Popule a parte Y\n");
    for(i=0;i<numeroPontosTabelados;i++){
        printf("y%d: ",i);
        scanf("%lf",&tabelaPontos[1][i]);
    }

    int selecaoFuncao;
    int grau=0;
    printf("selecione a funcao:\n");
    printf("1- newton\n");
    printf("2- newton gregory\n");
    printf("3- ajuste reta\n");
    printf("4- ajuste polinomio\n");
    printf("5- ajuste exponencial\n");
    printf("selecione a funcao:\n");
    scanf("%d",&selecaoFuncao);

    if(selecaoFuncao == 4){
        printf("informe o grau do polinomio\n");
        scanf("%d",&grau);      
    }
    double vetorCoeficientes[grau+1];

    double pontoDesejado;
    switch (selecaoFuncao){
    case 1:
        printf("Ponto desejado:\n");
        scanf("%lf",&pontoDesejado);

        Newton(numeroPontosTabelados,tabelaPontos,pontoDesejado);
        break;
    case 2:
        printf("Ponto desejado:\n");
        scanf("%lf",&pontoDesejado);

        NewtonGregory(numeroPontosTabelados,tabelaPontos,pontoDesejado);
        break;
    case 3:
        AjusteReta(numeroPontosTabelados,tabelaPontos,a0,a1,yajustado,coefDeterminacao);
        break;
    case 4:
        AjustePolinomial(numeroPontosTabelados,grau,tabelaPontos,vetorCoeficientes,yajustado,coefDeterminacao);
        break;
    case 5:
        AjusteExponencial(numeroPontosTabelados,tabelaPontos,a,b,yajustado,coefDeterminacao);
        break;
    
    default:
        break;
    }

    return (EXIT_SUCCESS);
}