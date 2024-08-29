# Needleman Wunsh MPI

<div style="justify-content: center"> 
  <img src="https://64.media.tumblr.com/072ce8ed7a1ae906be3d48edd4c5faff/6f9d4a3279d367fa-94/s1280x1920/b9c6823d25f72adccf59e81f905658456b283521.gif">
</div>

## Objetivo

O objetivo deste trabalho prático é desenvolver uma versão MPI da aplicação de reconhecimento de genoma fornecida pelo professor. **Não é permitida a utilização de outra aplicação de base, apenas a fornecida pelo professor.**

## Diretrizes

1. **Permitir que o usuário defina os tamanhos das sequências a serem comparadas, que poderão ser da ordem de milhares de bases.**
2. **Permitir que as sequências de bases possam ser geradas aleatoriamente, fornecidas em arquivos ou lidas a partir do teclado, todas as três opções.** Na geração aleatória, gera-se a primeira sequência aleatória e com base nesta, cria-se a segunda sequência. A segunda sequência deve ser retirada da primeira e deve sofrer um percentual de alterações, para ficarem com certo grau de similaridade. O percentual de alterações deve ser fornecido pelo usuário. No fornecimento por arquivo ou na leitura pelo teclado, as duas sequências devem ser fornecidas pelo usuário.
3. **As sequências devem ser compostas por caracteres que representam as bases (A, T, G e C).**
4. **A matriz de pesos deve ser fornecida pelo usuário via digitação.**
5. **Deve-se usar a mesma Função de Pontuação do material passado em aula, exceto pelo valor da penalidade do gap (d), que deve ser fornecido pelo usuário.**
6. **O programa deve mostrar o maior escore para o alinhamento global.** Além disso, o programa deverá mostrar até np possibilidades de alinhamento global e as sequências alinhadas pareadas uma com a outra.
7. **A matriz de scores deve ser gravada em arquivo texto, de forma tabulada e bem autoexplicativa.** As duas sequências de entrada devem aparecer nas linhas e colunas da matriz de scores de forma a permitir analisar e conferir a corretude da matriz.
8. **Deve-se entregar a solução MPI, além de um relatório técnico.** A versão MPI deve possuir as mesmas funcionalidades e características da versão sequencial.
9. **O relatório técnico deve informar os dados dos autores e explicar o funcionamento dos códigos, as restrições e limitações, os erros, as partes que não funcionam e como os códigos devem ser usados, com exemplos de uso.**
10. **Na solução MPI, o usuário deve informar o número de processos (np) quando da execução.**
11. **O paralelismo deve ser realizado em apenas uma etapa da aplicação: a construção/preenchimento da matriz de scores deve ser paralelizada entre np processos, de forma mais igualitária possível.** Tal distribuição deverá ser feita por linhas de forma que cada processo deve construir um conjunto de linhas. Cada processo deve trabalhar com um conjunto de linhas equidistantes, não consecutivas, para promover a paralelização. O traceback será único, identificando apenas um único alinhamento, executado pelo processo 0.
12. **Para permitir o paralelismo na construção da matriz de escores, na medida em que um processo vai gerando os escores de uma linha, ele vai transmitindo esta linha, bloco por bloco, para o próximo processo que precisa dela e para o processo 0 também.** Essa transmissão (por mensagem MPI) deve ser feita por blocos, para não congestionar muito o subsistema de troca de mensagens. **Os processos não devem esperar a finalização de uma linha para somente depois poder transmiti-la, senão a execução será sequencializada.**
13. **O tamanho do bloco deverá ser definido pelo usuário.** Como o processo 0 vai recebendo os blocos de todos os demais processos, ao final, ele terá a matriz de escores montada por completo e poderá assim executar o traceback. Nesse modelo de paralelização, o processo 0 não precisa trabalhar na paralelização, apenas no gerenciamento da execução MPI.
14. **A equipe de trabalho deve ser composta por no máximo 3 participantes.**

## Entrega

Cada equipe deverá criar um repositório no GitHub para este trabalho prático e fazer o push dos arquivos fonte em linguagem C padrão (.c) e do relatório técnico contendo as seguintes informações:

- **Identificação dos participantes**
- **Descrição dos principais módulos desenvolvidos**
- **Auto-avaliação do funcionamento** (elencar as partes que funcionam corretamente, as partes que não funcionam corretamente e sob quais circunstâncias, bem como as partes que não foram implementadas).

---

**Data Máxima de Entrega: 08/09/2024 - 23h55min**
