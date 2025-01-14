# Diferenças Finitas 2D
A seguir segue uma descrição da utilização do algoritmo, que implementa uma solução numérica para a equação de onda em 2 dimensões, a partir do método de diferenças finitas.

## Execução
O algoritmo implementado aplica o método de diferenças finitas, de forma que a cada vez que o número de iteração temporal atingir um múltiplo da variável *frame* o atual estado da grade será exportado para um arquivo binário, localizado na pasta *"output/multi_matrix"*.

Ademais, o algoritmo cria um arquivo de informações, que contém, além do tempo de execução, dados que serão utilizados por *"plot.py"* para realizar a produção do GIF de forma automática. 

## Makefile
O arquivo Makefile busca garantir a correta execução do programa em diferentes sistemas, para isso asseguramos que as devidas pastas serão criadas como pré-requisito para a execução do programa.

Além disso, ele conta com as seguintes sub-rotinas:
- `make clean`: Realiza a limpeza dos arquivos `.o`.
- `make run`: Compila e executa o arquivo, `algorithm.c`, em seguida realiza a chamada do arquivo `plot.py`, que faz os gráficos e os salva na pasta *"figures"*. Por fim, executa `make clean` para remover os arquivos objeto.
- `make graph`: Caso o arquivo *"2D\_wave\_equation.gif"* já exista, abre a imagem no reprodutor padrão do sistema. Caso contrário, executa `make run` e em seguida abre a imagem pra visualização.
