# 🧬 Needleman-Wunsch MPI

## 📝 Descrição
Este projeto implementa o algoritmo de Needleman-Wunsch para alinhamento global de sequências genômicas utilizando MPI (Message Passing Interface) para paralelização. O programa oferece uma interface interativa com diversas opções para manipulação e análise de sequências genômicas.

## ✨ Funcionalidades
- 📊 Leitura e exibição da matriz de pesos
- 🔢 Definição da penalidade de gap
- 🧬 Entrada de sequências genômicas:
  - ⌨️ Manual
  - 🎲 Geração aleatória
  - 📁 Leitura de arquivo
- 🔢 Geração e exibição da matriz de escores
- 🔍 Geração de alinhamento global:
  - 🥇 Primeiro maior
  - 🏁 Último maior
- 👀 Exibição do alinhamento global

## 🛠️ Requisitos
- 💻 Compilador C
- 🌐 Biblioteca MPI

## 🚀 Compilação
Para compilar o programa, use o seguinte comando no terminal:
```
mpicc MPI_NW.c -o MPI_NW
```

## 🏃‍♂️ Execução
Para executar o programa, use o seguinte comando:
```
mpirun -np <número de processos> MPI_NW <opções>
```
> 💡 Substitua `<número de processos>` pelo número desejado de processos.

## 🎮 Uso
Ao executar o programa, você verá um menu interativo com as seguintes opções:

1. 📥 Ler Matriz de Pesos
2. 📊 Mostrar Matriz de Pesos
3. 🔢 Ler Penalidade de Gap
4. 👁️ Mostrar Penalidade
5. 🧬 Definir Sequências Genômicas
6. 📜 Mostrar Sequências
7. 🔢 Gerar Matriz de Escores
8. 📊 Mostrar Matriz de Escores
9. 🔍 Gerar Alinhamento Global
10. 👀 Mostrar Alinhamento Global
11. 🚪 Sair

> 📌 Siga as instruções na tela para utilizar cada funcionalidade.