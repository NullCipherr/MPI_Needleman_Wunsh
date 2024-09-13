#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <limits.h>

#define A 0 // representa uma base Adenina
#define T 1 // representa uma base Timina
#define G 2 // representa uma base Guanina
#define C 3 // representa uma base Citosina
#define X 4 // representa um gap

#define sair 11

#define maxSeq 1000 // tamanho maximo de bases em uma sequencia genomica

void inicializaMatrizEscores();
void calculaElementos(int idSend, int idRec);
void verificaRecebimento(int idRec, int lin, int col);
void calculaElementoAtual(int lin, int col);
void enviaBloco(int idSend, int lin, int col, int *aux_send);
void enviaElementosFaltantes(int idSend, int lin, int aux_send);
void recebeBlocos();
void calculaMaioresElementos();

void traceBack(int tipo);
int performTraceback(int tbLin, int tbCol, int pos);
int getMaxIndex(int escores[3]);
int applyMovement(int maxIndex, int tbLin, int tbCol, int pos);
int insertRemainingGaps(int tbLin, int tbCol, int pos);
void reverseAlignmentSequences();

int np; // Numero de processos
int id; // Id do processo
int blocoS; // Bloco de sequencia

/* mapaBases mapeia indices em caracteres que representam as bases, sendo 0='A',
1='T', 2='G', 3='C' e 4='-' representando gap */

char mapaBases[5] = {'A', 'T', 'G', 'C', '-'};

int seqMaior[maxSeq] = {A, A, C, T, T, A},
    seqMenor[maxSeq] = {A, C, T, T, G, A};

int  alinhaGMaior[maxSeq],
     alinhaGMenor[maxSeq];

MPI_Status st;

int maxT;

int matrizEscores[maxSeq + 1][maxSeq + 1];

int tamSeqMaior = 6, /* tamanho da sequencia maior, inicializado como 6 */
    tamSeqMenor = 6, /* tamanho da sequencia menor, inicializado como 6 */
    tamAlinha,       /* tamanho do alinhamento global obtido */
    penalGap = 0,    /* penalidade de gap, a ser descontada no escore acumulado
                        quando um gap eh encontrado */
    grauMuta = 0,    /* porcentagem maxima de mutacao na geracao aleatoria da
                        sequencia menor, a qual eh copiada da maior e sofre algumas
                        trocas de bases */
    escoreDiag,      /* escore da diagonal anterior da matriz de escores */
    escoreLin,       /* escore da linha anterior da matriz de escores */
    escoreCol;       /* escore da coluna anterior da matriz de escores */

int matrizPesos[4][4] = {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1}};

int indRef = -1,                  // indice da sequencia maior a partir do qual extrai a sequencia
                                  // menor, no caso de geracao aleatoria
    nTrocas = -1,                 // quantidade de trocas na geracao automatica da sequencia menor,
                                  // a partir de um segmento da sequencia maior
    linPMaior, colPMaior, PMaior, // suporte para deteccao do primeiro maior escore
    linUMaior, colUMaior, UMaior; // suporte para deteccao do ultimo maior escore


// Função para imprimir mensagens de debug de forma profissional
void debug(const char *message) {
    printf("╔════════════════════════════════════════════════════════╗\n");
    printf("║ DEBUG: %-48s║\n", message);
    printf("╚════════════════════════════════════════════════════════╝\n");
}


// Função para mostrar o alinhamento global
void mostraAlinhamentoGlobal(void)
{   
    int i; // Variável auxiliar para o índice

    // Imprime o alinhamento obtido
    printf("\nAlinhamento Obtido - Tamanho = %d:\n", tamAlinha);
    fflush(stdout);

    printf("%c",mapaBases[alinhaGMaior[0]]); // Imprime a primeira base da sequência maior alinhada
    fflush(stdout);
    // Imprime a sequência maior alinhada
    for (i=1; i<tamAlinha; i++)
    {
        printf("%c",mapaBases[alinhaGMaior[i]]); // Imprime a sequência maior alinhada
        fflush(stdout);
    }
    printf("\n");
    fflush(stdout);

    printf("%c",mapaBases[alinhaGMenor[0]]); // Imprime a primeira base da sequência menor alinhada
    fflush(stdout);
    // Imprime a sequência menor alinhada
    for (i=1; i<tamAlinha; i++)
    {
        printf("%c",mapaBases[alinhaGMenor[i]]); // Imprime a sequência menor alinhada
        fflush(stdout);
    }
    printf("\n");
    fflush(stdout);
}


// Função para realizar o traceback
void traceBack(int tipo)
{
    int tbLin, tbCol, pos = 0;

    // Determina o tipo de alinhamento a ser gerado (Primeiro ou Último Maior)
    if (tipo == 1) 
    {
        debug("Geracao do Primeiro Maior Alinhamento Global");
        fflush(stdout);
        tbLin = linPMaior;
        tbCol = colPMaior;
    } 
    else 
    {
        debug("Geracao do Ultimo Maior Alinhamento Global");
        fflush(stdout);
        tbLin = linUMaior;
        tbCol = colUMaior;
    }

    // Realiza o traceback até que uma das sequências seja completamente percorrida
    pos = performTraceback(tbLin, tbCol, pos);

    // Define o tamanho do alinhamento
    tamAlinha = pos;

    // Inverte as sequências alinhadas para a ordem correta
    reverseAlignmentSequences();

    // Imprime mensagem de debug indicando a conclusão do alinhamento
    debug("Geracao do Alinhamento Global Concluida");
    fflush(stdout);
}


// Função para realizar o traceback 
int performTraceback(int tbLin, int tbCol, int pos)
{
    int movimentos[3][2] = {{-1, -1}, {0, -1}, {-1, 0}}; // Diagonal, esquerda, cima

    // Enquanto a linha e a coluna não forem 0, continua o traceback
    while (tbLin != 0 && tbCol != 0) 
    {
        int escores[3] = {
            matrizEscores[tbLin - 1][tbCol - 1] + matrizPesos[seqMenor[tbLin - 1]][seqMaior[tbCol - 1]], // Diagonal
            matrizEscores[tbLin][tbCol - 1] - penalGap, // Esquerda
            matrizEscores[tbLin - 1][tbCol] - penalGap  // Cima
        };

        int maxIndex = getMaxIndex(escores); // Obtém o índice do maior escore

        pos = applyMovement(maxIndex, tbLin, tbCol, pos); // Aplica o movimento

        tbLin += movimentos[maxIndex][0]; // Atualiza a linha
        tbCol += movimentos[maxIndex][1]; // Atualiza a coluna
    }

    pos = insertRemainingGaps(tbLin, tbCol, pos); // Insere os gaps restantes

    return pos; // Retorna a posição final
}


// Função para obter o índice do maior escore
int getMaxIndex(int escores[3])
{
    int maxIndex = 0; // Inicializa o índice do maior escore

    // Encontra o índice do maior escore
    for (int i = 1; i < 3; i++) 
    {
        // Se o escore atual é maior que o máximo, atualiza o índice do maior escore
        if (escores[i] > escores[maxIndex]) 
        {
            maxIndex = i; // Atualiza o índice do maior escore
        }
    }
    return maxIndex; // Retorna o índice do maior escore
}


// Função para aplicar o movimento
int applyMovement(int maxIndex, int tbLin, int tbCol, int pos)
{
    // Aplica o movimento de acordo com o índice do maior escore
    if (maxIndex == 0) 
    {
        // Se as bases são diferentes, insere gaps nas duas sequências
        if (seqMenor[tbLin - 1] != seqMaior[tbCol - 1]) 
        {
            printf("\nALERTA no TraceBack: Pos = %d Lin = %d e Col = %d\n", pos, tbLin, tbCol);
            fflush(stdout);
            alinhaGMenor[pos] = X;
            alinhaGMaior[pos] = seqMaior[tbCol - 1];
            pos++;
            alinhaGMenor[pos] = seqMenor[tbLin - 1];
            alinhaGMaior[pos] = X;
            pos++;
        } 
        else 
        {
            alinhaGMenor[pos] = seqMenor[tbLin - 1];
            alinhaGMaior[pos] = seqMaior[tbCol - 1];
            pos++;
        }
    } 
    else if (maxIndex == 1) // Esquerda
    {
        alinhaGMenor[pos] = X;
        alinhaGMaior[pos] = seqMaior[tbCol - 1];
        pos++;
    } 
    else // Cima
    {
        alinhaGMenor[pos] = seqMenor[tbLin - 1];
        alinhaGMaior[pos] = X;
        pos++;
    }
    return pos;
}


// Função para inserir os gaps restantes
int insertRemainingGaps(int tbLin, int tbCol, int pos)
{
    // Insere os gaps restantes da sequência menor
    while (tbLin > 0) 
    {
        alinhaGMenor[pos] = seqMenor[tbLin - 1]; // Insere a base da sequência menor
        alinhaGMaior[pos] = X; // Insere um gap na sequência maior
        tbLin--; // Decrementa o índice da linha
        pos++; // Incrementa a posição
    }

    // Insere os gaps restantes da sequência maior
    while (tbCol > 0) 
    {
        alinhaGMenor[pos] = X; // Insere um gap na sequência menor
        alinhaGMaior[pos] = seqMaior[tbCol - 1]; // Insere a base da sequência maior
        tbCol--; // Decrementa o índice da coluna
        pos++; // Incrementa a posição
    }

    return pos; // Retorna a posição final
}


// Função para inverter as sequências alinhadas
void reverseAlignmentSequences()
{
    int aux; // Variável auxiliar para armazenar o valor temporário
    
    // Inverte as sequências alinhadas
    for (int i = 0; i < tamAlinha / 2; i++) 
    {
        aux = alinhaGMenor[i]; // Armazena o valor atual
        alinhaGMenor[i] = alinhaGMenor[tamAlinha - i - 1]; // Inverte o valor
        alinhaGMenor[tamAlinha - i - 1] = aux; // Armazena o valor invertido

        aux = alinhaGMaior[i]; // Armazena o valor atual
        alinhaGMaior[i] = alinhaGMaior[tamAlinha - i - 1]; // Inverte o valor
        alinhaGMaior[tamAlinha - i - 1] = aux; // Armazena o valor invertido
    }
}


// Função para ler o tamanho da sequencia maior
void leTamMaior(void)
{
    debug("Leitura do Tamanho da Sequencia Maior");
    fflush(stdout);

    // Loop para garantir que o tamanho da sequencia maior seja válido
    do
    {
        printf("\nDigite 0 < valor < %d = ", maxSeq);
        fflush(stdout);
        scanf("%d", &tamSeqMaior);
    } while ((tamSeqMaior < 1) || (tamSeqMaior > maxSeq));
}


// Função para ler o tamanho da sequencia menor
void leTamMenor(void)
{
    debug("Leitura do Tamanho da Sequencia Menor");
    fflush(stdout);

    // Loop para garantir que o tamanho da sequencia menor seja válido
    do
    {
        printf("\nDigite 0 < valor <= %d = ", tamSeqMaior);
        fflush(stdout);
        scanf("%d", &tamSeqMenor);
    } while ((tamSeqMenor < 1) || (tamSeqMenor > tamSeqMaior));
}


// Função para ler a penalidade de gap
int lePenalidade(void)
{
    int penal; // Variável para armazenar a penalidade de gap

    debug("Leitura da Penalidade de Gap");
    fflush(stdout);

    // Loop para garantir que a penalidade de gap seja válida
    do
    {
        printf("\nDigite valor >= 0 = ");
        fflush(stdout);
        scanf("%d", &penal);
    } while (penal < 0);

    return penal;
}


// Função para ler a matriz de pesos
void leMatrizPesos()
{
    int i, j; // Variáveis auxiliares para os índices da matriz

    debug("Leitura da Matriz de Pesos");
    fflush(stdout);

    // Loop para ler a matriz de pesos
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            printf("Digite valor %c x %c = ", mapaBases[i], mapaBases[j]);
            fflush(stdout);
            scanf("%d", &(matrizPesos[i][j]));
        }
        printf("\n");
        fflush(stdout);
    }
}

// Função para mostrar a matriz de pesos
void mostraMatrizPesos(void)
{
    int i, j;

    debug("Imprimindo a Matriz de Pesos");
    fflush(stdout);
    printf("\n%4c%4c%4c%4c%4c\n", ' ', 'A', 'T', 'G', 'C');
    fflush(stdout);
    for (i = 0; i < 4; i++)
    {
        printf("%4c", mapaBases[i]);
        fflush(stdout);
        for (j = 0; j < 4; j++)
        {
            printf("%4d", matrizPesos[i][j]);
            fflush(stdout);
        }
        printf("\n");
        fflush(stdout);
    }
}


// Função para ler a porcentagem máxima de mutação
int leGrauMutacao(void)
{
    int prob;

    debug("Leitura da Porcentagem Maxima de Mutacao");
    fflush(stdout);

    // Loop para garantir que a porcentagem máxima de mutação seja válida
    do
    {
        printf("\nDigite 0 <= valor <= 100 = ");
        fflush(stdout);
        scanf("%d", &prob);
    } while ((prob < 0) || (prob > 100));

    return prob; // Retorna a porcentagem máxima de mutação
}


// Função para ler as sequências
void leSequencias(void)
{
    int i, erro;
    char seqMaiorAux[maxSeq], seqMenorAux[maxSeq];

    indRef = -1;
    nTrocas = -1;
    
    debug("Leitura das Sequências");
    fflush(stdout);

    /* lendo a sequencia maior */
    do
    {
        printf("\nPara a Sequencia Maior,");
        fflush(stdout);
        printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
        fflush(stdout);
        do
        {
            printf("\n> ");
            fflush(stdout);
            fgets(seqMaiorAux, maxSeq, stdin);
            tamSeqMaior = strlen(seqMaiorAux) - 1; /* remove o enter */
        } while (tamSeqMaior < 1);
        printf("\ntamSeqMaior = %d\n", tamSeqMaior);
        fflush(stdout);
        i = 0;
        erro = 0;
        do
        {
            switch (seqMaiorAux[i])
            {
            case 'A':
                seqMaior[i] = A;
                break;
            case 'T':
                seqMaior[i] = T;
                break;
            case 'G':
                seqMaior[i] = G;
                break;
            case 'C':
                seqMaior[i] = C;
                break;
            default:
                erro = 1; /* nao eh permitido qquer outro caractere */
            }
            i++;
        } while ((erro == 0) && (i < tamSeqMaior));
    } while (erro == 1);

    /* lendo a sequencia menor */
    do
    {
        printf("\nPara a Sequencia Menor, ");
        fflush(stdout);
        printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
        fflush(stdout);
        do
        {
            printf("\n> ");
            fflush(stdout);
            fgets(seqMenorAux, maxSeq, stdin);
            tamSeqMenor = strlen(seqMenorAux) - 1; /* remove o enter */
        } while ((tamSeqMenor < 1) || (tamSeqMenor > tamSeqMaior));
        
        printf("\ntamSeqMenor = %d\n", tamSeqMenor);
        fflush(stdout);

        i = 0;
        erro = 0;
        do
        {
            switch (seqMenorAux[i])
            {
            case 'A':
                seqMenor[i] = A;
                break;
            case 'T':
                seqMenor[i] = T;
                break;
            case 'G':
                seqMenor[i] = G;
                break;
            case 'C':
                seqMenor[i] = C;
                break;
            default:
                erro = 1;
            }
            i++;
        } while ((erro == 0) && (i < tamSeqMenor));
    } while (erro == 1);
}


// Função para gerar sequências aleatórias
void geraSequencias(void)
{
    int i, dif, probAux;
    char base;

    debug("Geracao Aleatoria das Sequencias");
    fflush(stdout);

    /* gerando a sequencia maior */
    for (i = 0; i < tamSeqMaior; i++)
    {
        base = rand() % 4; /* produz valores de 0 a 3 */
        seqMaior[i] = base;
    }

    dif = tamSeqMaior - tamSeqMenor; /* diferenca entre os tamanhos das sequencias */

    indRef = 0;
    if (dif > 0)
        indRef = rand() % dif; /* produz um indice aleatorio para indexar a sequencia maior,
                               para a partir dele extrair a primeira versao da sequencia
                               menor */

    /* gerando a sequencia menor a partir da maior. Copia trecho da sequencia
       maior, a partir de um indice aleatorio que nao ultrapasse os limites do
       vetor maior */
    for (i = 0; i < tamSeqMenor; i++)
        seqMenor[i] = seqMaior[indRef + i];

    /* causa mutacoes aleatorias na sequencia menor para gerar "gaps",
       sobre cada base, de acordo com o grau (porcentagem) informado.
       A mutacao causa a troca da base original por outra base aleatoria
       necessariamente diferente. Gera-se uma probabilidade aleatoria
       ateh 100 e se ela estiver dentro do grau informado, a mutacao
       ocorre na base, caso contrario, mantem a base copiada. */

    i = 0;
    nTrocas = 0;
    while ((i < tamSeqMenor) && (nTrocas < ((grauMuta * tamSeqMenor) / 100)))
    {
        probAux = rand() % 100 + 1;

        if (probAux <= grauMuta)
        {
            seqMenor[i] = (seqMenor[i] + (rand() % 3) + 1) % 4;
            nTrocas++;
        }
        i++;
    }

    printf("\nSequencias Geradas: Dif = %d, IndRef = %d, NTrocas = %d\n ", dif, indRef, nTrocas);
    fflush(stdout);
}


// Função para mostrar as sequências
void mostraSequencias(void)
{
    int i;

    debug("Imprimindo as Sequências    ");
    fflush(stdout);
    printf("\nSequencia Maior, Tam = %d\n", tamSeqMaior);
    fflush(stdout);

    for (i = 0; i < tamSeqMaior; i++)
    {
        printf("%c", mapaBases[seqMaior[i]]);
        fflush(stdout);
    }

    printf("\n");
    fflush(stdout);

    for (i = 0; i < tamSeqMaior; i++)
    {
        if (i != indRef)
        {
            printf(" ");
        }
        else
        {
            printf("^");
        }
    }
    printf("\nIndice de Referencia = %d\n", indRef);
    fflush(stdout);

    printf("\nSequencia Menor, Tam = %d\n", tamSeqMenor);
    fflush(stdout);

    for (i = 0; i < tamSeqMenor; i++)
    {
        printf("%c", mapaBases[seqMenor[i]]);
        fflush(stdout);
    }
    printf("\n");
    fflush(stdout);

    for (i = 0; i < tamSeqMenor; i++){
        if (seqMenor[i] != seqMaior[indRef + i])
        {
            printf("^");
            fflush(stdout);
        }
        else{
            printf(" ");
            fflush(stdout);
        }
    }
    printf("\nQuantidade de trocas = %d\n", nTrocas);
    fflush(stdout);
}


// Função para mostrar a matriz de escore
void mostraMatrizEscores()
{
    FILE *file;
    int i, lin, col;

    file = fopen("Matriz_Escores.txt", "w");

    if (file == NULL)
    {
        printf("ERRO ao criar o arquivo !!!\n");
        fflush(stdout);
        exit(1);
    }

    fprintf(file, "Matriz de escores Atual:\n");
    debug("Imprimindo a Matriz de Escores");
    fflush(stdout);

    fprintf(file, "%4c%4c", ' ', ' ');
    printf("%4c%4c", ' ', ' ');

    for (i = 0; i <= tamSeqMaior; i++)
    {
        fprintf(file, "%4d", i);
        printf("%4d", i);
        fflush(stdout);
    }

    fprintf(file, "\n");
    printf("\n");
    fflush(stdout);

    fprintf(file, "%4c%4c%4c", ' ', ' ', '-');
    printf("%4c%4c%4c", ' ', ' ', '-');
    fflush(stdout);

    for (i = 0; i < tamSeqMaior; i++)
    {
        fprintf(file, "%4c", mapaBases[(seqMaior[i])]);
        printf("%4c", mapaBases[(seqMaior[i])]);
        fflush(stdout);
    }

    fprintf(file, "\n");
    printf("\n");
    fflush(stdout);

    fprintf(file, "%4c%4c", '0', '-');
    printf("%4c%4c", '0', '-');
    fflush(stdout);

    for (col = 0; col <= tamSeqMaior; col++)
    {
        fprintf(file, "%4d", matrizEscores[0][col]);
        printf("%4d", matrizEscores[0][col]);
        fflush(stdout);
    }

    fprintf(file, "\n");
    printf("\n");
    fflush(stdout);

    for (lin = 1; lin <= tamSeqMenor; lin++)
    {
        fprintf(file, "%4d%4c", lin, mapaBases[(seqMenor[lin - 1])]);
        printf("%4d%4c", lin, mapaBases[(seqMenor[lin - 1])]);
        fflush(stdout);

        for (col = 0; col <= tamSeqMaior; col++)
        {
            fprintf(file, "%4d", matrizEscores[lin][col]);
            printf("%4d", matrizEscores[lin][col]);
            fflush(stdout);
        }
        fprintf(file, "\n");
        printf("\n");
        fflush(stdout);
    }

    fclose(file);
    debug("Matriz_Escores.txt salvo com sucesso");
    fflush(stdout);
}


// Função para gerar a matriz de escore
void geraMatrizEscores() 
{
    int idSend = id + 1; // id do processo que vai enviar os blocos
    int idRec = id - 1; // id do processo que vai receber os blocos

    if(idRec == 0) 
    {
        idRec = np - 1; // id do processo que vai receber os blocos
    }

    // Se o id do processo que vai enviar os blocos é maior que o número de processos, ele volta para 1
    if(idSend >= np)
    {
        idSend = 1; // id do processo que vai enviar os blocos
    }

    inicializaMatrizEscores(); // Inicializa a matriz de escore

    // Se o id do processo não é 0, ele calcula os elementos da matriz de escore
    if(id != 0)
    {
        calculaElementos(idSend, idRec);
    }

    // Se o id do processo é 0, ele recebe os blocos e calcula os maiores elementos da matriz de escore
    if(id == 0)
    {
        recebeBlocos(); // Recebe os blocos
        calculaMaioresElementos(); // Calcula os maiores elementos da matriz de escore
    }
}


// Função para inicializar a matriz de escore
void inicializaMatrizEscores()
{
    for(int i = 0; i <= tamSeqMenor; i++)
    {
        matrizEscores[i][0] =  -1 * (i * penalGap); // Inicializa a primeira coluna com o valor da penalidade de gap
    }

    for(int col = 0; col <= tamSeqMaior; col++ )
    {
        matrizEscores[0][col] = -1 * (col * penalGap); // Inicializa a primeira linha com o valor da penalidade de gap
    }
}


// Função para calcular os elementos da matriz de escore
void calculaElementos(int idSend, int idRec)
{
    // Calcula os elementos da matriz de escore
    for(int lin = id; lin <= tamSeqMenor; lin+= (np- 1))
    {
        int aux_send = 0; // aux_send é a coluna que falta enviar

        // Calcula os elementos da matriz de escore
        for(int col = 1; col <= tamSeqMaior; col++)
        {
            verificaRecebimento(idRec, lin, col); // Verifica se o elemento foi recebido
            calculaElementoAtual(lin, col); // Calcula o elemento atual

            // Envia o bloco se a coluna é múltiplo de blocoS
            if(col % blocoS == 0)
            {
                enviaBloco(idSend, lin, col, &aux_send);
            }
            
            // Envia os elementos faltantes se a coluna é maior que o tamanho da sequência maior
            if(col == tamSeqMaior && aux_send != tamSeqMaior)
            {
                enviaElementosFaltantes(idSend, lin, aux_send);
            }
        }
    }
}


// Função para verificar se o elemento foi recebido
void verificaRecebimento(int idRec, int lin, int col)
{
    if(matrizEscores[lin-1][col] == INT_MIN && col + blocoS <= tamSeqMaior) // Se o elemento foi recebido
    {
        MPI_Recv(&matrizEscores[lin-1][col], blocoS, MPI_INT, idRec, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Recebe o bloco
    }
    if(matrizEscores[lin-1][col] == INT_MIN && col + blocoS >= tamSeqMaior + 1) // Se o elemento foi recebido
    {
        MPI_Recv(&matrizEscores[lin-1][col], tamSeqMaior - col + 1, MPI_INT, idRec, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Recebe o bloco
    }
}


// Função para calcular o elemento atual
void calculaElementoAtual(int lin, int col)
{
    int peso = matrizPesos[(seqMenor[lin-1])][(seqMaior[col-1])]; // Calcula o peso da base
    int escoreDiag = matrizEscores[lin-1][col-1] + peso; // Calcula o escore da diagonal
    int escoreLin = matrizEscores[lin][col-1] - penalGap; // Calcula o escore da linha
    int escoreCol = matrizEscores[lin-1][col] - penalGap; // Calcula o escore da coluna

    if ((escoreDiag > escoreLin) && (escoreDiag > escoreCol)) // Se o escore da diagonal é maior que o escore da linha e da coluna
    {
        matrizEscores[lin][col] = escoreDiag; // Atualiza o escore da célula
    }
    else if (escoreLin > escoreCol) // Se o escore da linha é maior que o escore da coluna
    {
        matrizEscores[lin][col] = escoreLin; // Atualiza o escore da célula
    }
    else // Se o escore da coluna é maior que o escore da linha
    {
        matrizEscores[lin][col] = escoreCol; // Atualiza o escore da célula
    }
}


// Função para enviar os blocos
void enviaBloco(int idSend, int lin, int col, int *aux_send)
{
    // Envia o bloco para o processo idSend
    if(lin != tamSeqMenor)
    {
        MPI_Send(&matrizEscores[lin][col-blocoS + 1], blocoS, MPI_INT, idSend, 0, MPI_COMM_WORLD); // Envia o bloco para o processo idSend
    }

    // Envia o tamanho do bloco, a linha e a coluna para o processo 0
    int colSend = col-blocoS + 1; // colSend é a coluna que falta enviar
    MPI_Send(&blocoS, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // Envia o tamanho do bloco para o processo 0
    MPI_Send(&lin, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // Envia a linha para o processo 0
    MPI_Send(&colSend, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // Envia a coluna para o processo 0
    MPI_Send(&matrizEscores[lin][colSend], blocoS, MPI_INT, 0, 0, MPI_COMM_WORLD); // Envia o bloco para o processo 0
    *aux_send = col; // Atualiza o aux_send com o valor da coluna
}


// Função para enviar os elementos faltantes
void enviaElementosFaltantes(int idSend, int lin, int aux_send)
{
    int sendSize = tamSeqMaior - aux_send; // sendSize é o tamanho do bloco que falta enviar
    int colSend = aux_send + 1; // colSend é a coluna que falta enviar
    if(lin != tamSeqMenor) // Se a linha não for a última linha
    {
        MPI_Send(&matrizEscores[lin][aux_send + 1], sendSize, MPI_INT, idSend, 0, MPI_COMM_WORLD); // Envia o bloco que falta enviar
    }
    MPI_Send(&sendSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // Envia o tamanho do bloco que falta enviar 
    MPI_Send(&lin, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // Envia a linha que falta enviar
    MPI_Send(&colSend, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // Envia a coluna que falta enviar
    MPI_Send(&matrizEscores[lin][colSend], sendSize, MPI_INT, 0, 0, MPI_COMM_WORLD); // Envia o bloco que falta enviar
}


// Função para receber os blocos
void recebeBlocos()
{
    int lin, col; // lin e col são a linha e a coluna que receberão o bloco

    // Recebe os blocos que faltam enviar
    for(int i = 1; i < np; i++)
    {
        for(int l = i; l <= tamSeqMenor; l += (np - 1))
        {
            for(int c = 1; c <= tamSeqMaior; c += blocoS)
            {
                int recSize = 0; // recSize é o tamanho do bloco que falta receber
                MPI_Recv(&recSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Recebe o tamanho do bloco que falta receber
                MPI_Recv(&lin, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Recebe a linha que falta receber
                MPI_Recv(&col, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Recebe a coluna que falta receber
                MPI_Recv(&matrizEscores[lin][col], recSize, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Recebe o bloco que falta receber
            }
        }
    }
}


// Função para calcular os maiores elementos da matriz de escore
void calculaMaioresElementos()
{
    linPMaior = 1; // linPMaior é a linha do primeiro maior escore
    colPMaior = 1; // colPMaior é a coluna do primeiro maior escore 
    PMaior = matrizEscores[1][1]; // PMaior é o primeiro maior escore

    linUMaior = 1; // linUMaior é a linha do último maior escore
    colUMaior = 1; // colUMaior é a coluna do último maior escore
    UMaior = matrizEscores[1][1]; // UMaior é o último maior escore

    // Calcula os maiores elementos da matriz de escore
    for (int lin = 1; lin <= tamSeqMenor; lin++)
    {
        for (int col = 1; col <= tamSeqMaior; col++)
        {
            if (PMaior < matrizEscores[lin][col])
            {
                linPMaior = lin; // linPMaior é a linha do primeiro maior escore
                colPMaior = col; // colPMaior é a coluna do primeiro maior escore
                PMaior = matrizEscores[lin][col]; // PMaior é o primeiro maior escore
            }
            if (UMaior <= matrizEscores[lin][col])
            {
                linUMaior = lin; // linUMaior é a linha do último maior escore
                colUMaior = col; // colUMaior é a coluna do último maior escore
                UMaior = matrizEscores[lin][col]; // UMaior é o último maior escore
            }
        }
    }

    debug("Matriz de escores Gerada.");
    fflush(stdout);
    printf("\nPrimeiro Maior escore = %d na celula [%d,%d]", PMaior, linPMaior, colPMaior);
    fflush(stdout);
    printf("\nUltimo Maior escore = %d na celula [%d,%d]", UMaior, linUMaior, colUMaior);
    fflush(stdout);
}


// Função para mostrar o menu de opções
int menuOpcao(void)
{
    int op;
    char enter;

    do
    {
        printf("\n");
        printf("╔════════════════════════════════════════════════════════╗\n");
        printf("║                  Needleman Wunsh MPI                   ║\n");
        printf("╠════════════════════════════════════════════════════════╣\n");
        printf("║  01. Ler Matriz de Pesos                               ║\n");
        printf("║  02. Mostrar Matriz de Pesos                           ║\n");
        printf("║  03. Ler Penalidade de Gap                             ║\n");
        printf("║  04. Mostrar Penalidade                                ║\n");
        printf("║  05. Definir Sequências Genômicas                      ║\n");
        printf("║  06. Mostrar Sequências                                ║\n");
        printf("║  07. Gerar Matriz de Escores                           ║\n");
        printf("║  08. Mostrar Matriz de Escores                         ║\n");
        printf("║  09. Gerar Alinhamento Global                          ║\n");
        printf("║  10. Mostrar Alinhamento Global                        ║\n");
        printf("║  11. Sair                                              ║\n");
        printf("╚════════════════════════════════════════════════════════╝\n");
        printf("╔════════════════════════════════════════════════════════╗\n");
        printf("║ Digite a opção => ");
        fflush(stdout);
        scanf("%d", &op);
        scanf("%c", &enter);
        fflush(stdout);
        printf("╚════════════════════════════════════════════════════════╝\n");
        fflush(stdout);
    } while ((op < 1) || (op > sair));

    return op;
}


// Função para ler as sequências de um arquivo
void leArquivo() 
{
    // Declaração de variáveis para armazenar o nome do arquivo e as sequências
    char nome_arquivo[100];
    char seq_maior[maxSeq], seq_menor[maxSeq];
    int erro = 0;

    // Solicita ao usuário o caminho do arquivo
    printf("Digite o caminho do arquivo: \n");
    fflush(stdout);
    fgets(nome_arquivo, sizeof(nome_arquivo), stdin);
    nome_arquivo[strcspn(nome_arquivo, "\n")] = 0;  // Remove newline

    // Abre o arquivo para leitura
    FILE *file = fopen(nome_arquivo, "r");

    // Verifica se o arquivo foi aberto corretamente
    if (file == NULL) 
    {
        perror("Erro ao abrir o arquivo !!!");
        fflush(stdout);
        exit(1);
    }

    // Lê 2 linhas do arquivo, correspondentes às 2 sequências
    if (fgets(seq_maior, sizeof(seq_maior), file) == NULL || fgets(seq_menor, sizeof(seq_menor), file) == NULL) {
        perror("Erro ao ler as sequências do arquivo !!!");
        fflush(stdout);
        fclose(file);
        exit(1);
    }

    // Remove o caractere de nova linha das sequências lidas
    seq_maior[strcspn(seq_maior, "\n")] = 0;
    seq_menor[strcspn(seq_menor, "\n")] = 0;

    // Define os tamanhos das sequências
    tamSeqMaior = strlen(seq_maior);
    tamSeqMenor = strlen(seq_menor);

    // Fecha o arquivo
    fclose(file);

    // Converte caracteres da seq_maior para os valores correspondentes
    for (int i = 0; i < tamSeqMaior; i++) 
    {
        switch (seq_maior[i]) 
        {
            case 'A':
                seqMaior[i] = A;
                break;
            case 'T':
                seqMaior[i] = T;
                break;
            case 'G':
                seqMaior[i] = G;
                break;
            case 'C':
                seqMaior[i] = C;
                break;
            default:
                erro = 1;
                debug("Erro: A sequencia maior contem caracteres invalidos\n");
                fflush(stdout);
                break;
        }
        if (erro) break; // Interrompe a conversão se houver erro
    }

    // Converte caracteres da seq_menor para os valores correspondentes
    if (!erro) {
        for (int i = 0; i < tamSeqMenor; i++) 
        {
            switch (seq_menor[i]) {
            case 'A':
                seqMenor[i] = A;
                break;
            case 'T':
                seqMenor[i] = T;
                break;
            case 'G':
                seqMenor[i] = G;
                break;
            case 'C':
                seqMenor[i] = C;
                break;
            default:
                erro = 1;
                debug("Erro: A sequencia menor contem caracteres invalidos\n");
                fflush(stdout);
                break;
            }
            if (erro) break; // Interrompe a conversão se houver erro
        }
    }

    // Se não houver erros, exibe as sequências lidas com sucesso
    if (!erro) 
    {
        debug("A sequencia maior e menor foram lidas com sucesso!\n");
        fflush(stdout);
        mostraSequencias();
    }
}


// Função para inicializar a matriz de escores
void initializeScoreMatrix()
{
    // Inicializa a matriz de escores com o valor INT_MIN
    for(int i = 0; i < maxSeq; i++)
    {
        for(int j = 0; j < maxSeq; j++)
        {
            matrizEscores[i][j] = INT_MIN; // INT_MIN é o menor valor inteiro
        }
    }

    // Se o processo é o mestre (id == 0)
    if(id == 0)
    {
        // Solicita ao usuário o tamanho do bloco
        printf("\nDigite o tamanho do bloco: ");
        fflush(stdout);
        scanf("%d", &blocoS);

        // Envia o tamanho do bloco, tamanho das sequências e penalidade de gap para todos os outros processos
        for(int i = 1; i < np; i++)
        {
            MPI_Send(&blocoS, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&tamSeqMaior, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&tamSeqMenor, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&penalGap, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        // Gera a matriz de escores
        geraMatrizEscores();
    } 
    else 
    {
        // Se o processo não é o mestre, recebe os dados do processo mestre
        MPI_Recv(&blocoS, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&tamSeqMaior, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&tamSeqMenor, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&penalGap, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Gera a matriz de escores
        geraMatrizEscores();
    }
}

// Função para mostrar os créditos do programa
void mostraCreditos() 
{
    printf("\n\n");
    printf("╔════════════════════════════════════════════════════════╗\n");
    printf("║                                                        ║\n");
    printf("║            Programa Needleman-Wunsch MPI               ║\n");
    printf("║                                                        ║\n");
    printf("║            Desenvolvido por:                           ║\n");
    printf("║            Andrei Costa     - RA107975                 ║\n");
    printf("║            João Casagrande  - RA112684                 ║\n");
    printf("║            Kananda Caroline - RA116382                 ║\n");
    printf("║                                                        ║\n");
    printf("║            Agradecimentos especiais a:                 ║\n");
    printf("║            Prof. Ronaldo Goncalves                     ║\n");
    printf("║            Universidade Estadual de Maringá            ║\n");
    printf("║                                                        ║\n");
    printf("╚════════════════════════════════════════════════════════╝\n");
    printf("\n\n");
}


// Função para tratar a opção selecionada pelo usuário
void trataOpcao(int op)
{
    int resp; // Variável para armazenar a resposta do usuário
    char enter; // Variável para armazenar o enter

    // Trata a opção selecionada
    switch (op)
    {
        case 1: // Lê a matriz de pesos
            debug("Lendo Matriz de Pesos");
            fflush(stdout);
            leMatrizPesos();
            break;
        case 2: // Mostra a matriz de pesos
            debug("Mostrando Matriz de Pesos");
            mostraMatrizPesos();
            break;
        case 3: // Lê a penalidade de gap
            debug("Lendo Penalidade de Gap");
            penalGap = lePenalidade();
            break;
        case 4: // Mostra a penalidade de gap
            debug("Mostrando Penalidade de Gap");
            printf("\nPenalidade = %d", penalGap);
            break;
        case 5: // Define as sequências genômicas
            debug("Definindo Sequências Genômicas  ");
            do {
                printf("╔════════════════════════════════════════════════════════╗\n");
                printf("║               Definição de Sequências Genômicas        ║\n");
                printf("╠════════════════════════════════════════════════════════╣\n");
                printf("║  <1> MANUAL                                            ║\n");
                printf("║  <2> ALEATÓRIA                                         ║\n");
                printf("║  <3> ARQUIVO                                           ║\n");
                printf("║  <4> VOLTAR                                            ║\n");
                printf("╚════════════════════════════════════════════════════════╝\n");
                printf("╔════════════════════════════════════════════════════════╗\n");
                printf("║ Digite a opção => ");
                fflush(stdout);
                if (scanf("%d", &resp) != 1) {
                    fprintf(stderr, "Erro: Entrada inválida\n");
                    break;
                }
                scanf("%c", &enter); // Consome o enter
                printf("╚════════════════════════════════════════════════════════╝\n");
                fflush(stdout);
                
                if (resp == 1)
                {
                    leSequencias(); // Lê as sequências
                }
                else if (resp == 2) // Gera sequências aleatórias
                {
                    leTamMaior(); // Lê o tamanho da sequência maior
                    leTamMenor(); // Lê o tamanho da sequência menor
                    grauMuta = leGrauMutacao(); // Lê o grau de mutação
                    geraSequencias(); // Gera as sequências
                }
                else if (resp == 3)
                {
                    leArquivo(); // Lê as sequências de um arquivo
                }
                else if (resp == 4)
                {
                    debug("Retornando ao menu principal");
                    return; // Volta ao menu principal
                }
                else
                {
                    fprintf(stderr, "Erro: Opção inválida\n");
                }
            } while (resp < 1 || resp > 4);
            break;

        case 6:
            debug("Mostrando Sequências  ");
            fflush(stdout);
            mostraSequencias();
            break;
        case 7:
            debug("Gerando Matriz de Escores");
            fflush(stdout);
            initializeScoreMatrix();
            break;
        case 8:
            debug("Mostrando Matriz de Escores");
            fflush(stdout);
            mostraMatrizEscores();
            break;
        case 9:
            debug("Gerando Alinhamento Global");
            fflush(stdout);

            printf("╔════════════════════════════════════════════════════════╗\n");
            printf("║ Deseja: <1> Primeiro Maior ou <2> Último Maior?        ║\n");

            printf("╚════════════════════════════════════════════════════════╝\n");

            printf("╔════════════════════════════════════════════════════════╗\n");

            printf("║ Digite a opção => ");
            fflush(stdout);

            if (scanf("%d", &resp) != 1) 
            {
                fprintf(stderr, "Erro: Entrada inválida\n");
                break;
            }
            scanf("%c", &enter); // Consome o enter

            printf("╚════════════════════════════════════════════════════════╝\n");
            fflush(stdout);

            traceBack(resp);
            break;
        case 10:
            debug("Mostrando Alinhamento Global");
            mostraAlinhamentoGlobal();
            break;
        case 11:
            debug("Saindo do programa...");
            mostraCreditos();
            break;
        default:
            debug("Erro: Opção inválida");
            break;
    }
}

// Função principal do programa
int main(int argc, char *argv[])
{
    // Variáveis locais
    int opcao = 0; // Variável para armazenar a opção selecionada pelo usuário

    // Inicializa a semente do gerador de números aleatórios
    srand(time(NULL));

    // Inicializa o ambiente MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id); // Obtém o rank do processo
    MPI_Comm_size(MPI_COMM_WORLD, &np); // Obtém o número total de processos

    // Loop principal do programa
    do 
    {
        if (id == 0) 
        {
            // Exibe o menu e obtém a opção do usuário
            opcao = menuOpcao();

            // Envia a opção selecionada para todos os outros processos
            for (int i = 1; i < np; i++) 
            {
                MPI_Send(&opcao, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }

            // Trata a opção selecionada
            trataOpcao(opcao);
        } 
        else 
        {
            // Processos escravos (id != 0)
            MPI_Recv(&opcao, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Se a opção for 7 (Gerar Matriz de Escores), prepara o processo
            if (opcao == 7) 
            {
                initializeScoreMatrix();
            }
        }
    } while (opcao != sair); // Continua até que a opção de sair seja selecionada

    // Finaliza o ambiente MPI
    MPI_Finalize();

    return 0;
}
